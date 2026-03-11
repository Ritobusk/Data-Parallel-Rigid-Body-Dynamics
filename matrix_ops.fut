-- Taken from https://futhark-lang.org/examples/matrix-multiplication.html
def matmul [n][m][p] 'a
           (add: a -> a -> a) (mul: a -> a -> a) (zero: a)
           (A: [n][m]a) (B: [m][p]a) : [n][p]a =
  map (\A_row ->
         map (\B_col ->
                reduce add zero (map2 mul A_row B_col))
             (transpose B))
      A
def matmul_f64 = matmul (+) (*) (0f64)

def mat_mul_vec [n] [m] 'a 
                (add: a -> a -> a) (mul: a -> a -> a) (zero: a)
                (A: [n][m]a) (v: [m]a) : [n]a =
  map (\A_row ->
          reduce add zero (map2 mul A_row v)
      ) A
def mat_mul_vec_f64 = mat_mul_vec (+) (*) (0f64)

def matadd [n][m] 'a
           (add : a -> a -> a) 
           (A: [n][m]a) (B: [n][m]a) : [n][m]a =
    tabulate_2d n m (\r c -> add A[r][c] B[r][c])
def matadd_f64 = matadd (f64.+) 

def scal_mul_mat [n] [m] (s : f64) (A : [n][m]f64) : [n][m]f64 =
    map (\r -> map (\x -> x * s) r) A 

def diagonal [a] (diag : [a]f64) : [a][a]f64 =
  tabulate_2d a a (\r c -> if r == c then diag[r] else 0f64)


def vec_mul_vec [n] (v1: [n]f64) (v2: [n]f64) : f64 =
  reduce (+) 0f64 (map2 (*) v1 v2)

-- Taken from https://futhark-lang.org/examples/swap.html
def swap 't (i: i64) (j: i64) (A: *[]t) =
  let tmp = copy A[j]
  let A[j] = copy A[i]
  let A[i] = tmp
  in A


-- Taken from https://futhark-lang.org/student-projects/kristian-kasper-peter-project.pdf
-- Work: O(n))
-- Span: O(log(n))
def argmax (arr: []f64) =
    reduce_comm (\(a, i) (b, j) ->
    if a < b
    then (b, j)
    else if b < a then (a, i)
    else if j < i then (b, j)
    else (a, i)
  ) (0, 0) (zip arr (indices arr))


-- Work: O(min(m, n) · m · n))
  -- Span: O(min(m, n) · log(m))
def gauss_jordan [m][n] (A:[m][n]f64) =
  loop A = copy A for i < i64.min m n do
    -- Find largest pivot
    let p = A[i:,i] |> map f64.abs |> argmax |> (.1) |> (+i)
    let A = if p != i then swap i p A else A
    let irow = map (/A[i,i]) A[i]
    -- Eliminate entries above and below the pivot
    in tabulate m (\j ->
        let scale = A[j,i]
        in map2 (\x y ->
          if j != i then y - scale * x else x
          ) irow A[j]
     )
  
-- let gauss_solveAB [m][n] (A:[m][m]f32) (B:[m][n]f32) : [m][n]f32 =
--    let AB = gauss_jordan (hStack A B)
--    in AB[:m, m:] :> [m][n]f32
--   
-- let gauss_solveAb [m] (A:[m][m]f32) (b:[m]f32) =
--    unflatten m 1 b |> gauss_solveAB A |> flatten_to m
