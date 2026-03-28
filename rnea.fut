import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"

module T = vtree

def mkt 'a [n] (ps:[n]i64) (ds:[n]a) : [n]{parent:i64,data:a} =
    map2 (\p d -> {parent=p,data=d}) ps ds
def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

-- Inspiration taken from https://royfeatherstone.org/spatial/v2/sourceText/ID.txt
-- As far as I understand q, qd and qdd are scalars. This might only be the case for joints with only 1-DOF
--  q, qd, qdd and tau are column vectors of length model.NB containing the joint position, velocity, acceleration and force variables, respectively.
def rnea [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd  S 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let (vs, as) = loop (vs', as') = (replicate n (replicate 6 0f64), replicate n (replicate 6 0f64)) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i]  ( (-1) `scal_mul_vec_f64` gravity  )) 
                                   (qdd[i] `scal_mul_vec_f64` S[i])
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = map2 (+) (Xup[i] `mat_mul_vec_f64`  vs'[parent])  vJ[i]

        let as' = as' with [i] = map2 (+) ( Xup[i] `mat_mul_vec_f64` as'[parent])   
                                           ( qdd[i]`scal_mul_vec_f64` S[i])
                                |> vecadd_f64 ( (crm vs'[i]) `mat_mul_vec_f64` vJ[i])
        in (vs', as')

  let fBs = map (\i -> 
              (mat_mul_vec_f64 Is[i] as[i]) `vecadd_f64` (((crf vs[i]) `matmul_f64` Is[i]) `mat_mul_vec_f64` vs[i])
              ) (iota n) 

  let (tau, _) = loop (tau', fs') = (replicate n 0f64, fBs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = vecmul S[idx]  fs'[idx] 
    in 
      if idx > 0 then
        let fs'' = fs' with [parent] = map2 (+) fs'[parent] ((transpose Xup[idx]) `mat_mul_vec_f64` fs'[idx])
        in (tau', fs'')
      else (tau', fs')
  in trace tau


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let vs = loop vs' = (copy vJ) for i < (n-1) do
      let parent = p[i+1]
      in vs' with [i+1] = map2 (+) (Xup[i+1] `mat_mul_vec_f64` vs'[parent]) 
                                    vs'[i+1]

  let as_tmp = map2 (vecadd_f64)
                        (map (\i -> qdd[i]  `scal_mul_vec_f64` S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity)
                          else  (crm vs[i]) `mat_mul_vec_f64` vJ[i]) (iota n))

  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
      let parent = p[i+1]
      in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) 
                                     as'[i+1]

  let fBs = map (\i -> 
              (Is[i] `mat_mul_vec_f64` as[i]) `vecadd_f64` (((crf vs[i]) `matmul_f64` Is[i]) `mat_mul_vec_f64` vs[i])
              ) (iota n) 

  let fJs = loop fJs' = (fBs) for i < n -1 do
      let idx = n - (i+1)
      let parent = p[idx]
      let fJs'' = fJs' with [parent] = map2 (+) (copy fJs'[parent])
                                    <| mat_mul_vec_f64 (transpose Xup[idx]) (copy fJs'[idx])
      in fJs''

  in map2 (\s f -> s `vecmul`  f) S fJs  


-- Same as above implementation but now with vtrees
def rnea'' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) =
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let Cs = zip Xup vJ

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_M ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let vs = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  -- let vtree_vs = T.mk_preorder <| mkt p Cs
  -- let vs2 = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let vs = map (.1) vs

  -- as_tmp = S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
  let as_tmp = map2 (\S_qdd v_cross_S_qd -> map2 (+) S_qdd v_cross_S_qd)
                        (map (\i -> map (\s -> s * qdd[i]) S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then replicate 6 0f64
                          else mat_mul_vec_f64 (crm vs[i]) vJ[i]) (iota n))
  let as_tmp = as_tmp with [0] = map2 (+) as_tmp[0] (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))

  let Cs = zip Xup as_tmp

  let vtree_as = T.lprp <| mkt2 lp rp Cs
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
  let as = map (.1) as -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  -- Ide 1: Får de 2 nederste niveauer til at være korrekt, men man får ikke "+ f_Bi" termet 
  --        med til resten af niveauerne.
  -- let Cs = map2 (\m v -> (transpose m) `mat_mul_vec_f64` v) Xup fBs 
  --
  -- let inv_op (ci : [6]f64) : [6]f64 =
  --    scal_mul_vec_f64 (-1) (ci)
  --
  -- let operator (si : [6]f64) (ci : [6]f64) : [6]f64 =
  --   si `vecadd_f64` ci
  -- let vtree_fs = T.lprp <| mkt2 lp rp Cs
  -- let fJs2 = T.ileaffix_sc operator inv_op (replicate 6 0f64) vtree_fs

  -- let Cs = zip (map (transpose) Xup) fBs
  --
  -- let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
  --   let inv_cia = identity 6
  --   let inv_cib = scal_mul_vec_f64 (-1) (ci.1)
  --   in (inv_cia, inv_cib)
  --
  -- let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
  --   (ci.0,  (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
  --
  --
  -- let vtree_fs = T.lprp <| mkt2 lp rp Cs
  -- let fJs2 = T.ileaffix_sc operator inv_op (identity 6, replicate 6 0f64) vtree_fs
  -- let fJs2 = map (.1) fJs2

  let vtree_transformation = T.lprp <| mkt2 lp rp Xup
  let from_root_M = T.irootfix (matmul_f64) (XBtoA_from_XAtoB_M) (identity 6) vtree_transformation

  let to_root_F   = map (transpose) from_root_M
  let from_root_F = map XBtoA_MtoF from_root_M

  let fBs_root = map2 (\i_to_root fbi -> i_to_root `mat_mul_vec_f64` fbi) to_root_F fBs

  let inv_op (ci : [6]f64) : [6]f64 =
    scal_mul_vec_f64 (-1) (ci)

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 inv_op (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (\X_to_joint fji -> X_to_joint `mat_mul_vec_f64` fji) from_root_F fJs_root 

  in map2 (\s f -> vecmul s f) S fJs  


-- rnea with vtrees that also take an array of external forces
--  these external forces are expressed in root coordinates
-- To get an understanding of how this affects the total force see eq. 5.20 on page 95 in Roy Featherstones book
def rnea_vtree_with_f_ext [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (f_ext : [n][6]f64)
             (lp : [n]i64) (rp : [n]i64) =
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let Cs = zip Xup vJ

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_M ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let vs = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let vs = map (.1) vs

  let as_tmp = map2 (\S_qdd v_cross_S_qd -> map2 (+) S_qdd v_cross_S_qd)
                        (map (\i -> map (\s -> s * qdd[i]) S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then replicate 6 0f64
                          else mat_mul_vec_f64 (crm vs[i]) vJ[i]) (iota n))
  let as_tmp = as_tmp with [0] = map2 (+) as_tmp[0] (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))

  let Cs = zip Xup as_tmp

  let vtree_as = T.lprp <| mkt2 lp rp Cs
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
  let as = map (.1) as -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  let vtree_transformation = T.lprp <| mkt2 lp rp Xup
  let from_root_M = T.irootfix matmul_f64 XBtoA_from_XAtoB_M (identity 6) vtree_transformation
  let to_root_F   = map (transpose) from_root_M
  let from_root_F = map XBtoA_MtoF from_root_M

  -- Here you add the external forces
  --  Since you transform the body forces to root coordinates the external force can just be added 
  let fBs_root = map3 (\i_to_root fbi f_xi -> 
            (i_to_root `mat_mul_vec_f64` fbi)  `vecadd_f64` (scal_mul_vec_f64 (-1) f_xi)
              ) to_root_F fBs f_ext

  let inv_op (ci : [6]f64) : [6]f64 =
    scal_mul_vec_f64 (-1) (ci)

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 inv_op (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (\X_to_joint fji -> X_to_joint `mat_mul_vec_f64` fji) from_root_F fJs_root 

  in map2 (\s f -> vecmul s f) S fJs  


def main = 
  let lp = [0, 1, 5, 2]
  let rp = [7, 4, 6, 3]
  let (_, p, js, _, Is, Xtrees) = autoTree 4 2 1 1
  in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 1] [0f64, 2, 1, 3] [0f64, 3, 0,  3] lp rp
  -- let lp = [0, 1, 7, 2, 4, 8]
  -- let rp = [11,  6, 10, 3, 5, 9]
  -- let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  -- in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3] lp rp
  -- let (_, p, js, _, Is, Xtrees) = autoTree 100 1 1 1
  -- in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] (replicate 100 (1f64)) (replicate 100 (1f64)) (replicate 100 (1f64))
  -- in rnea' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3]
