
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

def skew (v : [3]f64) : [3][3]f64 =
  [[0,     -v[2], v[1]],
   [v[2],  0,     -v[0]],
   [-v[1], v[0],  0]]

-- Spacial translation from A to B
-- where r is the vector from A to B
def xlt (r : [3]f64) : [6][6]f64 =
    [[1,    0,    0,    0,0,0],
     [0,    1,    0,    0,0,0],
     [0,    0,    1,    0,0,0],
     [0,    r[2], -r[1],1,0,0],
     [-r[2],0,    r[0], 0,1,0],
     [r[1], -r[0],0,    0,0,1]]

-- Rotation of theta radians about the X-axis
def rotx (theta : f64) : [6][6]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[1,0, 0,0,0, 0],
       [0,c, s,0,0, 0],
       [0,-s,c,0,0, 0],
       [0,0, 0,1,0, 0],
       [0,0, 0,0,c, s],
       [0,0, 0,0,-s,c]]

-- Make rigid body inertia matrix
def mcI (m : f64) (CoM : [3]f64) (I : [3][3]f64) : [6][6]f64 =
    let C = skew CoM
    let C' = transpose C
    -- 'quodrants' of the 6x6 matrix
    let q2 = matadd_f64 I <| scal_mul_mat m <| matmul_f64 C C'
    let q1 = scal_mul_mat m C
    let q3 = scal_mul_mat m C'
    let q4 = diagonal [m,m,m]
    in tabulate_2d 6 6 
      (\r c -> if r < 3 then
                  if c < 3 then --q2
                    q2[r][c]
                  else          --q1
                    q1[r][c-3]
               else
                  if c < 3 then --q3
                    q3[r-3][c]
                  else          --q4
                    q4[r-3][c-3]
                )


type jointT = #revolute | #prismatic | #helical

-- Inspiration is taken from: https://royfeatherstone.org/spatial/v2/sourceText/autoTree.txt
-- Creates a kinematic tree 
def autoTree (nb : i64) (bf : f64) (skew : f64) (taper : f64) =
    let ids = trace <| iota nb
    let joint_types = replicate nb (#revolute) : [nb]jointT
    let parents = map (\i -> i64.f64 <| (f64.floor ( (((f64.i64 i) + 1.0) - 2.0 + (f64.ceil bf) ) / bf )) - 1.0 ) ids

    let lengths = map (\i -> taper ** (f64.i64 i)) ids 
    let masses = trace <| map (\i ->  taper ** (3 * (f64.i64 i))) ids 
    let CoMs   = map (\l -> [0.5 * l, 0, 0]) lengths
    let Icms   = map (\i ->
                  let m = masses[i]
                  let l = lengths[i]
                  let d = diagonal [0.0025,1.015/12,1.015/12]  -- each body is a cylinder
                  in scal_mul_mat (m * l) d) ids

    let Is = trace <| map (\i -> mcI masses[i] CoMs[i] Icms[i]) ids 

    let Xtree       = trace <| map (\i -> if i == 0 then xlt [0,0,0] -- identity
                                 else matmul_f64 (rotx skew)  (xlt [lengths[parents[i]], 0, 0]) -- Not sure if this just puts all the children in the excact same place or not...
                          ) ids

    in (ids, joint_types, parents, lengths, Is, Xtree)

def main =
  autoTree 2 1 1 1
