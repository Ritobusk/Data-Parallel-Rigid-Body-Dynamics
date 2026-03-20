import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"

module T = vtree

def mkt 'a [n] (ps:[n]i64) (ds:[n]a) : [n]{parent:i64,data:a} =
    map2 (\p d -> {parent=p,data=d}) ps ds


-- autoTree creates a tree where the first node is connected a fixed base.
-- All other nodes have bf children (the second parameter of autoTree) if the tree_size
-- allows it.

entry mk_preordertest (tree_size : i64)  =
  let children = 2
  let (_, p, _, _, _, _) = autoTree tree_size children 1 1
  let vtree =trace <|  T.mk_preorder <| mkt p (replicate tree_size 1)
  in p

def main =
  let t1 = trace <| mk_preordertest 4
  let t2 = trace <| mk_preordertest 6
  in t1
  -- Manually calculated lp and rp for tree of size 4
  --let lp = [0, 1, 5, 2]
  --let rp = [7, 4, 6, 3]
  -- Manually calculated lp and rp for tree of size 6
  -- let lp = [0, 1, 7, 2, 4, 8]
  -- let rp = [11,  6, 10, 3, 5, 9]
