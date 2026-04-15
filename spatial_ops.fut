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

-- Transform a motion Plucker transform to a force Plucker transfrom
def XBtoA_MtoF (XAtoB : [6][6]f64) : [6][6]f64 =
  tabulate_2d 6 6 
    (\r c -> if r < 3 then
                if c < 3 then --q2
                  XAtoB[r][c]
                else          --q1
                  XAtoB[r+3][c-3]
             else
                if c < 3 then --q3
                  0f64
                else          --q4
                  XAtoB[r][c]
              )
 
-- Motion
-- Given a transformation from A to B get to inverse, i.e. 
-- the transformation from B to A. This can be done by composing
-- a transformation of the inverse rotation of A to B and the
-- 'inverse' vector of r (vector locating B in A). 
def XBtoA_from_XAtoB_M (XAtoB : [6][6]f64) : [6][6]f64 =
  let inv_E  = transpose XAtoB[:3, :3]

  let rotated_rx = sized 3 (XAtoB[3:6, :3])
  let org_rx = matmul_f64 inv_E rotated_rx -- this is: -rx as seen from A
  let inv_r  = matmul_f64 (scal_mul_mat_f64 (-1) org_rx) (inv_E)

  in tabulate_2d 6 6 
    (\r c -> if r < 3 then
                if c < 3 then --q2
                  inv_E[r][c]
                else          --q1
                  0f64
             else
                if c < 3 then --q3
                  inv_r[r-3][c]
                else          --q4
                  inv_E[r-3][c-3]
              )

-- Force version of XBtoA_from_XAtoB_M 
def XBtoA_from_XAtoB_F (XAtoB : [6][6]f64) : [6][6]f64 =
  let inv_E  = sized 3 <| transpose XAtoB[:3, :3]

  let rotated_rx = sized 3 (XAtoB[:3, 3:6])
  let org_rx =  matmul_f64 inv_E rotated_rx :> [3][3]f64 -- this is: -rx as seen from A
  let inv_r  = matmul_f64 (scal_mul_mat_f64 (-1) org_rx) (inv_E)

  in tabulate_2d 6 6 
    (\r c -> if r < 3 then
                if c < 3 then --q2
                  inv_E[r][c]
                else          --q1
                  inv_r[r][c-3]
             else
                if c < 3 then --q3
                  0f64
                else          --q4
                  inv_E[r-3][c-3]
              )


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
  crm v |> transpose |> scal_mul_mat_f64 (-1f64) 

-- Make rigid body inertia matrix
def mcI (m : f64) (CoM : [3]f64) (I : [3][3]f64) : [6][6]f64 =
    let C = skew CoM
    let C' = transpose C
    -- 'quodrants' of the 6x6 matrix
    let q2 = matadd_f64 I <| scal_mul_mat_f64 m <| matmul_f64 C C'
    let q1 = scal_mul_mat_f64 m C
    let q3 = scal_mul_mat_f64 m C'
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
