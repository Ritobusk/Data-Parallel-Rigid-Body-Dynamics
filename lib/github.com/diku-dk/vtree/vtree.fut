-- vtrees - data-parallel implementation of trees based on Euler-tours (Tarjan
-- and Vishkin) and Blelloch's insights (see "Guy Blelloch. Vector Models for Data-Parallel
-- Computing, MIT Press, 1990" (https://www.cs.cmu.edu/~guyb/papers/Ble90.pdf)

module type vtree = {
  type t 'a [n]

  val mk_preorder 'a [n] : [n]{parent: i64, data: a} -> t a [n]

  -- preorder node numbering
  val lprp 'a [n] : {lp: [n]i64, rp: [n]i64, data: [n]a} -> t a [n]

  -- preorder node numbering

  val map 'a 'b [n] : (a -> b) -> t a [n] -> t b [n]

  val rootfix 'a [n] :
    (op: a -> a -> a)
    -> (inv: a -> a)
    -> (ne: a)
    -> t a [n] -> [n]a

  val irootfix 'a [n] :
    (op: a -> a -> a)
    -> (inv: a -> a)
    -> (ne: a)
    -> t a [n] -> [n]a

  val leaffix 'a [n] :
    (op: a -> a -> a)
    -> (inv: a -> a)
    -> (ne: a)
    -> t a [n] -> [n]a

  val ileaffix 'a [n] :
    (op: a -> a -> a)
    -> (inv: a -> a)
    -> (ne: a)
    -> t a [n] -> [n]a

  val depth 'a [n] : t a [n] -> [n]i64
}

-- [mk_preorder a] creates a vtree from the preorder specification `a` of a tree
-- where each node is specified with data and a parent pointer (index into the
-- preorder of the nodes in the tree).
--
-- [lprp a] creates a vtree from a preorder specification `a` of a tree where
-- the `lp` and `rp` arrays specify the indices of a left- and right-parenthesis
-- print of the tree.
--
-- [map f t] maps a function `f` over the nodes in the tree.
--
-- [rootfix f inv ne t] computes, for each node `n` in the tree with the path
-- `n0->n1->...->nm->n` from the root `n0` to `n`, the values `f n0 (f n1
-- (...nm...))` (i.e., excluding the data for node `n`).
--
-- [rootfixi f inv ne t] computes, for each node `n` in the tree with the path
-- `n0->n2->...->n` from the root `n0` to `n`, the values `f n0 (f n1
-- (...n...))` (i.e., including the data for node `n`).
--
-- [leaffix f inv ne t] computes, for each node `n` in the tree with descendants
-- `[n1,n2,...,nm]` the accumulated value `f n1 (f n2 (...nm...))` (i.e.,
-- excluding the data for node `n`).
--
-- [leaffixi f inv ne t] computes, for each node `n` in the tree with
-- descendants `[n1,n2,...,nm]` the accumulated value `f n (f n1 (f n2
-- (...nm...)))` (i.e., including the data for node `n`).
--
-- For all `Xfix` operations (i.e., rootfix, rootfixi, leaffix, leaffixi), `f`
-- is assumed to be associative, the value `ne` is supposed to be the neutral
-- element for `f`, and `inv` is assumed to be the inverse associated with `f`
-- such that if `c = f a b` then `a = f c (inv b)` and `b = f c (inv a)`. All
-- `Xfix` operations has Work O(`m`) and Depth O(1), assuming `f` has constant
-- work and depth complexities.
--
-- [depth t] returns, for each node `n`, the depth of the subtree rooted at `n`
-- in `t`.

module vtree : vtree = {
  -- A vtree is represented by its left-parenthesis array and its
  -- right-parenthesis array.
  type t 'a [n] = {lp: [n]i64, rp: [n]i64, data: [n]a}

  type t0 'a [n] = [n]{parent: i64, data: a}

  -- A tree consists of n nodes (vertices, named 0..n-1) and n-1 edges. An edge
  -- p->x is uniquely identified by x (the node pointed to). There is no edge 0.
  --
  -- For constructing a vtree, we first construct an Euler tour of the tree and then
  -- extract the lp and rp arrays.

  def wyllie_scan_step [n] 'a
                       (op: a -> a -> a)
                       (values: [n]a)
                       (parents: [n]i64) =
    let f i =
      if parents[i] == -1
      then (values[i], parents[i])
      else (values[i] `op` values[parents[i]], parents[parents[i]])
    in unzip (tabulate n f)

  def wyllie_scan [n] 'a
                  (op: a -> a -> a)
                  (values: [n]a)
                  (parents: [n]i64) =
    let (values, _) =
      loop (values, parents) for _i < 64 - i64.clz n do
        wyllie_scan_step op values parents
    in values

  def exscan f ne xs =
    map2 (\i x -> if i == 0 then ne else x)
         (indices xs)
         (rotate (-1) (scan f ne xs))

  def size (h: i64) : i64 =
    (1 << h) - 1

  def mk_tree [n] 't (op: t -> t -> t) (ne: t) (arr: [n]t) =
    let temp = i64.num_bits - i64.clz n
    let h = i64.i32 <| if i64.popc n == 1 then temp else temp + 1
    let tree_size = size h
    let offset = size (h - 1)
    let offsets = iota n |> map (+ offset)
    let tree = scatter (replicate tree_size ne) offsets arr
    let arr = copy tree[offset:]
    let (tree, _, _) =
      loop (tree, arr, level) = (tree, arr, h - 2)
      while level >= 0 do
        let new_size = length arr / 2
        let new_arr =
          tabulate new_size (\i -> arr[2 * i] `op` arr[2 * i + 1])
        let offset = size level
        let offsets = iota new_size |> map (+ offset)
        let new_tree = scatter tree offsets new_arr
        in (new_tree, new_arr, level - 1)
    in tree

  def find_next [n] 't
                (op: t -> t -> bool)
                (tree: [n]t)
                (idx: i64) : i64 =
    let sibling i = i - i64.bool (i % 2 == 0) + i64.bool (i % 2 == 1)
    let parent i = (i - 1) / 2
    let is_right i = i % 2 == 0
    let h = i64.i32 <| i64.num_bits - i64.clz n
    let offset = size (h - 1)
    let start = offset + idx
    let v = tree[start]
    let ascent i = i != 0 && (is_right i || !(tree[sibling i] `op` v))
    let descent i = 2 * i + 2 - i64.bool (tree[2 * i + 1] `op` v)
    let index = iterate_while ascent parent start
    in if index != 0
       then iterate_while (< offset) descent (sibling index) - offset
       else -1

  def mk_preorder 'a [n] (ns: t0 a [n]) : t a [n] =
    let data = map (.data) ns
    let parents = map (.parent) ns
    -- Compute depths of each node in tree.
    let ds = wyllie_scan (+) (tabulate n (i64.bool <-< (!= 0))) parents
    -- Every node can be seen as edge on the euler tour in the
    -- downwards direction, if the depth changes then it meant that
    -- some subpath is taken up the tree.  These are the missing edges
    -- in the euler path. The upwards going edges and downwards going
    -- edges can be seen as parenthesis so these are the number of
    -- missing right parenthesis after an given left parenthesis.
    let missing =
      map3 (\i d d' -> if i != n - 1 && d >= d' then d - d' + 1 else 0)
           (indices ds)
           ds
           (rotate 1 ds)
    -- Adjust left parenthesis indices to account for the number of
    -- right parenthesis that are needed to be added.
    let lp = map2 (+) (indices ds) (exscan (+) 0 missing)
    -- Scatter the left parenthesis to their new position.
    let parens = scatter (replicate (2 * n) (-1)) lp (rep 1)
    -- Compute the depth of every parenthesis.
    let depths =
      parens
      |> scan (+) 0
      |> map2 (\p d -> d - i64.bool (p == 1)) parens
    -- Construct a prefix tree of minima.
    let min_tree = mk_tree i64.min i64.highest depths
    -- Find the next smaller or equal element.
    let rp = map (find_next (<=) min_tree) lp
    in {lp, rp, data}

  def lprp 'a [n] (x: t a [n]) : t a [n] = x

  def rootfix 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) ({lp, rp, data}: t a [n]) : [n]a =
    let lp = trace lp
    let rp = trace rp
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = exscan op ne R
    in map (\i -> S[i]) lp

  def irootfix 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (t: t a [n]) : [n]a =
    map2 op (rootfix op inv ne t) t.data

  def ileaffix 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) ({lp, rp, data}: t a [n]) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let S = exscan op ne L
    let Rv = map (\i -> S[i]) rp
    let Lv = map (\i -> inv (S[i])) lp
    in map2 op Rv Lv

  def leaffix 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (t: t a [n]) : [n]a =
    map2 op (ileaffix op inv ne t) (map inv t.data)

  def map 'a 'b [n] (f: a -> b) ({lp, rp, data}: t a [n]) : t b [n] =
    {lp, rp, data = map f data}

  def depth 'a [n] (t: t a [n]) : [n]i64 =
    let t' = map (\_ -> 1) t
    in rootfix (i64.+) i64.neg 0 t'
}
