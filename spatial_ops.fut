import "matrix_ops"

type jointT = #Rx | #Ry | #Rz | #Px | #Py | #Pz | #helical f64

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

-- Rotation of theta radians about the Y-axis
def roty (theta : f64) : [6][6]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[c,0,-s,0,0, 0],
       [0,1, 0,0,0, 0],
       [s,0, c,0,0, 0],
       [0,0, 0,c,0,-s],
       [0,0, 0,0,1, 0],
       [0,0, 0,s,0,c]]

-- Rotation of theta radians about the Y-axis
def rotz (theta : f64) : [6][6]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[c, s, 0,0, 0, 0],
       [-s,c, 0,0, 0, 0],
       [0, 0, 1,0, 0, 0],
       [0, 0, 0,c, s,0],
       [0, 0, 0,-s,c, 0],
       [0, 0, 0,0, 0,1]]

--  crm  spatial/planar cross-product operator (force).
--  crm(v)*f is the cross product of the motion vectors v and force f.
def crm (v : [6]f64) : [6][6]f64 =
  [ [ 0f64, -v[2],  v[1],   0   ,  0   ,  0    ], 
    [ v[2],  0   , -v[0],   0   ,  0   ,  0    ],
    [-v[1],  v[0],  0   ,   0   ,  0   ,  0    ],
    [ 0   , -v[5],  v[4],   0   , -v[2],  v[1] ],
    [ v[5],  0   , -v[3],   v[2],  0   , -v[0] ],
    [-v[4],  v[3],  0   ,  -v[1],  v[0],  0 ]]

--  crf  spatial/planar cross-product operator (motion).
--  crm(v)*m is the cross product of the motion vectors v and m.
def crf (v : [6]f64) : [6][6]f64 =
  crm v |> transpose |> scal_mul_mat (-1f64) 

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

-- Reverse of mcI
-- function  [m,c,I] = rbi_to_mcI( rbi )
--
-- if all(size(rbi)==[6 6])		% spatial
--
--   m = rbi(6,6);
--   mC = rbi(1:3,4:6);
--   c = skew(mC)/m;
--   I = rbi(1:3,1:3) - mC*mC'/m;


-- For all revolute joints, and also the helical joint, q is an angle in radians.  For all prismatic joints, q is a length in metres.  (If you know what you are doing, then you can choose a different length unit; but you must be sure you are using a consistent set of physical units overall.) 
def jcalc (jtyp : jointT) (q : f64) : ([6][6]f64, [6]f64) =
  match jtyp
  case  #Rx -> (rotx q, [1f64,0,0,0,0,0])
  case  #Ry -> (roty q, [0,1f64,0,0,0,0])
  case  #Rz -> (rotz q, [0,0,1f64,0,0,0])
  case  #Px -> (xlt [q,0,0], [0,0,0,1f64,0,0])
  case  #Py -> (xlt [0,q,0], [0,0,0,0,1f64,0])
  case  #Pz -> (xlt [0,0,q], [0,0,0,0,0,1f64])
  case  #helical pitch -> (replicate 6 (replicate 6 (0.0f64)), [0,0,1,0,0,pitch])
