import "matrix_ops"

type jointT = #Rx | #Ry | #Rz | #Px | #Py | #Pz | #helical f64
type mv = {w : [3]f64, v_O : [3]f64}
type fv = {n_O : [3]f64, f : [3]f64}
type X_Compact = {rot: [3][3]f64, r : [3]f64}
type I_Compact = {m: f64, h : [3]f64, I : [6]f64}

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
def rotx3d (theta : f64) : [3][3]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[1,0, 0],
       [0,c, s],
       [0,-s,c]]

-- Rotation of theta radians about the X-axis
def roty3d (theta : f64) : [3][3]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[c,0, -s],
       [0,1, 0],
       [s,0,c]]

-- Rotation of theta radians about the X-axis
def rotz3d (theta : f64) : [3][3]f64 =
    let s = f64.sin theta
    let c = f64.cos theta
    in 
      [[c, s, 0],
       [-s,c, 0],
       [0, 0, 1]]

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

-- Transform a motion Plucker transform to a force Plucker transfrom
def XBtoA_FtoM (XAtoB : [6][6]f64) : [6][6]f64 =
  tabulate_2d 6 6 
    (\r c -> if r < 3 then
                if c < 3 then --q2
                  XAtoB[r][c]
                else          --q1
                  0f64
             else
                if c < 3 then --q3
                  XAtoB[r-3][c+3]
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

def cross3d (v1 : [3]f64) (v2: [3]f64) : [3]f64 =
  [v1[1]*v2[2] - v1[2]*v2[1],
   v1[2]*v2[0] - v1[0]*v2[2],
   v1[0]*v2[1] - v1[1]*v2[0]]

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
def to_I_Compact (I : [6][6]f64) : I_Compact =
  let m = I[5,5]
  let mC = I[:3, 3:6]
  let c = [-mC[1,2], mC[0,2], -mC[0,1]]
  let Ic = I[:3, :3] --`matsub_f64` (scal_mul_mat_f64 (1/m) (mC `matmul_f64` (transpose mC)))
  let ltIc = lower_triangle_3d Ic 
  in {m = m, h = c, I = ltIc}


def IC_mul_mv (IC : I_Compact) (m : mv) : fv =
  let I = lt_unfold IC.I
  let n_O = (I `mat_mul_vec_f64` m.w) `vecadd_f64` (IC.h `cross3d` m.v_O)
  let f = ((IC.m) `scal_mul_vec_f64` m.v_O) `vecsub_f64` (IC.h `cross3d` m.w)
  in {n_O = n_O, f = f}

def fv_to_6d (f' : fv) : [6]f64 =
  sized 6 <| f'.n_O ++ f'.f

def d6_to_fv (f' : [6]f64) : fv =
  {n_O = sized 3 f'[0:3], f = sized 3 f'[3:6]}

-- constant times force vector
def scal_mul_fv (s : f64) (fv : fv) : fv =
  let n_O   = scal_mul_vec_f64 s  fv.n_O  
  let f     = scal_mul_vec_f64 s  fv.f  
  in {n_O = sized 3 n_O, f = sized 3 f}

-- add force vectors
def fv_add (fv1 : fv) (fv2 : fv) : fv =
  let n_O   = fv1.n_O `vecadd_f64`  fv2.n_O  
  let f     = fv1.f   `vecadd_f64`  fv2.f  
  in {n_O = n_O, f = f}

def mv_to_6d (m : mv) : [6]f64 =
  sized 6 <| m.w ++ m.v_O

def d6_to_mv (m : [6]f64) : mv =
  {w = sized 3 m[0:3], v_O = sized 3 m[3:6]}

-- constant times motion vector
def scal_mul_mv (s : f64) (mv : mv) : mv =
  let w   = scal_mul_vec_f64 s  mv.w  
  let v_O = scal_mul_vec_f64 s  mv.v_O  
  in {w = w, v_O = v_O}

-- add motion vector
def mv_add (mv1 : mv) (mv2 : mv) : mv =
  let w = mv1.w `vecadd_f64` mv2.w  
  let v_O = mv1.v_O `vecadd_f64` mv2.v_O 
  in {w = w, v_O = v_O}

def mv_cross_mv (mv1 : mv) (mv2 : mv) : mv =
  let w = mv1.w `cross3d` mv2.w  
  let v_O = (mv1.w `cross3d` mv2.v_O) `vecadd_f64` (mv1.v_O `cross3d` mv2.w) 
  in {w = w, v_O = v_O}

def mv_cross_fv (mv : mv) (fv : fv) : fv =
  let n_O = (mv.w `cross3d` fv.n_O) `vecadd_f64` (mv.v_O `cross3d` fv.f)
  let f = mv.w `cross3d` fv.f
  in {n_O = n_O, f = f}

def scalar_prod (m : mv) (f : fv) : f64 =
  (m.w `vecmul` f.n_O) + (m.v_O `vecmul` f.f)

-- From 6x6 Plücker transform matrix to compact representation
def to_X_Compact (X : [6][6]f64) : X_Compact =
  let E = X[:3, :3]
  let rotated_rx = sized 3 (X[3:6, :3])
  let org_rx = matmul_f64 (transpose E) rotated_rx -- this is: -rx as seen from A
  let r = [org_rx[1,2], -org_rx[0,2], org_rx[0,1]]
  in {rot = E, r = r}

-- From 6x6 Plücker transform matrix to compact representation
def from_XC_to_aXb (X : X_Compact) : [6][6]f64 =
  let E = X.rot
  let r = X.r
  let q3 = scal_mul_mat_f64 (-1) (matmul_f64 (E) (skew r) )
  in tabulate_2d 6 6 
    (\r c -> if r < 3 then
                if c < 3 then --q2
                  E[r][c]
                else          --q1
                  0
             else
                if c < 3 then --q3
                  q3[r-3][c]
                else          --q4
                  E[r-3][c-3]
              )

def X_inv (X : X_Compact) : X_Compact =
  let r = scal_mul_vec_f64 (-1) (X.rot `mat_mul_vec_f64` X.r)
  let E = transpose X.rot
  in {rot = E, r = r}

-- aXb Transform a motion vector
def Xm (X : X_Compact ) (m : mv) : mv =
  let w = X.rot `mat_mul_vec_f64` m.w
  let v_O = X.rot `mat_mul_vec_f64` (m.v_O `vecsub_f64` (X.r `cross3d` m.w) )
  in {w = w, v_O = v_O}

-- bXa Transform a motion vector
def Xm_inv (X : X_Compact ) (m : mv) : mv =
  let rot = transpose X.rot
  let w   = rot `mat_mul_vec_f64` m.w
  let v_O = (rot `mat_mul_vec_f64` m.v_O) `vecadd_f64` (X.r `cross3d` (w) )
  in {w = w, v_O = v_O}

-- aXb* Transform a force vector
def Xf (X : X_Compact ) (f : fv) : fv =
  let n_O = X.rot `mat_mul_vec_f64` (f.n_O `vecsub_f64` (X.r `cross3d` f.f) )
  let f   = X.rot `mat_mul_vec_f64` f.f
  in {n_O = n_O, f = f}

-- bXa* Transform a force vector
def Xf_inv (X : X_Compact ) (f : fv) : fv =
  let rot' = transpose X.rot
  let f'   = rot' `mat_mul_vec_f64` f.f
  let n_O = (rot' `mat_mul_vec_f64` f.n_O) `vecadd_f64` (X.r `cross3d` f') 
  in {n_O = n_O, f = f'}

def transform_XX (X1 : X_Compact ) (X2 : X_Compact) : X_Compact =
  let E = X1.rot `matmul_f64` X2.rot
  let r = X2.r `vecadd_f64` ((transpose X2.rot) `mat_mul_vec_f64` X1.r)
  in {rot = E, r = r}

def transform_XX_rev (X1 : X_Compact ) (X2 : X_Compact) : X_Compact = transform_XX X2 X1

def transform_identity : X_Compact =
  let E = identity 3
  let r  = [0,0,0f64]
  in {rot = E, r = r}

def transform_inv (X : X_Compact ) : X_Compact =
  let Et = transpose X.rot
  let r  = X.rot `mat_mul_vec_f64` (scal_mul_vec_f64 (-1) X.r)
  in {rot = Et, r = r}

-- For all revolute joints, and also the helical joint, q is an angle in radians.  For all prismatic joints, q is a length in metres.  (If you know what you are doing, then you can choose a different length unit; but you must be sure you are using a consistent set of physical units overall.) 
def jcalc (jtyp : jointT) (q : f64) : ([6][6]f64, [6]f64) =
  match jtyp
  case  #Rx -> (rotx q, [1f64,0,0,0,0,0])
  case  #Ry -> (roty q, [0,1f64,0,0,0,0])
  case  #Rz -> (rotz q, [0,0,1f64,0,0,0])
  case  #Px -> (xlt [q,0,0], [0,0,0,1f64,0,0])
  case  #Py -> (xlt [0,q,0], [0,0,0,0,1f64,0])
  case  #Pz -> (xlt [0,0,q], [0,0,0,0,0,1f64])
  case  #helical pitch -> 
      (rotz q `matmul_f64` xlt [0, 0, q * pitch], [0,0,1,0,0,pitch])

-- Only revolute joins are correct.
def jcalcC (jtyp : jointT) (q : f64) : (X_Compact, mv) =
  match jtyp
  case  #Rx -> ({rot = rotx3d q, r = [0,0,0f64]}, {w = [1f64,0,0], v_O = [0,0,0]})
  case  #Ry -> ({rot = roty3d q, r = [0,0,0f64]}, {w = [0,1f64,0], v_O = [0,0,0]})
  case  #Rz -> ({rot = rotz3d q, r = [0,0,0f64]}, {w = [0,0,1f64], v_O = [0,0,0]})
  case  #Px -> ({rot = rotz3d q, r = [0,0,0f64]}, {w = [0,0,1f64], v_O = [0,0,0]})
  case  #Py -> ({rot = rotz3d q, r = [0,0,0f64]}, {w = [0,0,1f64], v_O = [0,0,0]})
  case  #Pz -> ({rot = rotz3d q, r = [0,0,0f64]}, {w = [0,0,1f64], v_O = [0,0,0]})
  case  #helical pitch ->  ({rot = rotz3d q, r = [0,0,0f64]}, {w = [0,0,1f64], v_O = [0,0,0]})
