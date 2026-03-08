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

def matadd [n][m] 'a
           (add : a -> a -> a) 
           (A: [n][m]a) (B: [n][m]a) : [n][m]a =
    tabulate_2d n m (\r c -> add A[r][c] B[r][c])
def matadd_f64 = matadd (f64.+) 

def scal_mul_mat [n] [m] (s : f64) (A : [n][m]f64) : [n][m]f64 =
    map (\r -> map (\x -> x * s) r) A 

def diagonal [a] (diag : [a]f64) : [a][a]f64 =
  tabulate_2d a a (\r c -> if r == c then diag[r] else 0f64)

