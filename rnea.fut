import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"
import "scan_variations"

module T = vtree

def mkt 'a [n] (ps:[n]i64) (ds:[n]a) : [n]{parent:i64,data:a} =
    map2 (\p d -> {parent=p,data=d}) ps ds
def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

-- Inspiration taken from https://royfeatherstone.org/spatial/v2/sourceText/ID.txt
--  q, qd, qdd and tau are column vectors of length N containing the joint 
--   position, velocity, acceleration and force variables, respectively.
-- The parameters from the model that are used are the parent array, the joint types,
--  the inertia matrices and the Xtree transformations.
def rnea [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) : [n]f64 = 
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

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let (tau, _) = loop (tau', fs') = (replicate n 0f64, fBs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = vecmul S[idx]  fs'[idx] 
    in 
      if idx > 0 then
        let fs'' = fs' with [parent] = map2 (+) fs'[parent] ((transpose Xup[idx]) `mat_mul_vec_f64` fs'[idx])
        in (tau', fs'')
      else (tau', fs')
  in tau


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) 
             : [n]f64 = 
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let vs = loop vs' = (copy vJ) for i < (n-1) do
      let parent = p[i+1]
      in vs' with [i+1] = map2 (+) (Xup[i+1] `mat_mul_vec_f64` vs'[parent]) 
                                    vs'[i+1]

  let S_qdd  = tabulate n (\i -> scal_mul_vec_f64 qdd[i] S[i])
  let v_cross_S_qd = tabulate n (\i -> if i == 0 then (Xup[0] `mat_mul_vec_f64` (map (\x -> -1 * x) gravity))
                                                 else (crm vs[i]) `mat_mul_vec_f64` vJ[i])
  let as_tmp = map2 (vecadd_f64) S_qdd v_cross_S_qd 
  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
      let parent = p[i+1]
      in as' with [i+1] = map2 (+) (Xup[i+1] `mat_mul_vec_f64` as'[parent]) 
                                     as'[i+1]

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let fJs = loop fJs' = (fBs) for i < n -1 do
      let idx = n - (i+1)
      let parent = p[idx]
      let fJs'' = fJs' with [parent] = map2 (+) (copy fJs'[parent])
                                    <| mat_mul_vec_f64 (transpose Xup[idx]) (copy fJs'[idx])
      in fJs''

  in map2 (vecmul) S fJs  


-- Same as above implementation but now with vtrees
def rnea'' [n] (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : [n]f64 =
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64) XJ Xtree

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

  let S_qdd  = tabulate n (\i -> scal_mul_vec_f64 qdd[i] S[i])
  let v_cross_S_qd = tabulate n (\i -> if i == 0 then (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))
                                                 else mat_mul_vec_f64 (crm vs[i]) vJ[i])
  let as_tmp = map2 (vecadd_f64) S_qdd v_cross_S_qd 

  let Cs = zip Xup as_tmp

  let vtree_as = T.lprp <| mkt2 lp rp Cs
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
  let as = map (.1) as -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  let from_body_to_root_M = rootfix2 matmul_f64 XBtoA_from_XAtoB_M (identity 6) lp rp (map XBtoA_from_XAtoB_M Xup)

  let to_root_F   = map XBtoA_MtoF from_body_to_root_M 
  let from_root_F = map transpose from_body_to_root_M 

  let fBs_root = map2 (\X_to_root fbi -> X_to_root `mat_mul_vec_f64` fbi) to_root_F fBs

  let inv_op (ci : [6]f64) : [6]f64 =
    scal_mul_vec_f64 (-1) (ci)

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 inv_op (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (\X_to_joint fji -> X_to_joint `mat_mul_vec_f64` fji) from_root_F fJs_root 

  in map2 (vecmul) S fJs  

-- Same as above rnea'' but using different scans for vtrees
def rnea_vtree_optimized [n] (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : [n]f64 =
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64)       XJ Xtree

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_M ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) 
    = (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp (zip Xup vJ)
  let vs = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
          |> map (.1)

  let S_qdd  = tabulate n (\i -> scal_mul_vec_f64 qdd[i] S[i])
  let v_cross_S_qd = tabulate n (\i -> if i == 0 then (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))
                                                 else mat_mul_vec_f64 (crm vs[i]) vJ[i])
  let as_tmp = map2 (vecadd_f64) S_qdd v_cross_S_qd 

  let vtree_as = T.lprp <| mkt2 lp rp (zip Xup as_tmp)
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
          |> map (.1) 

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let vtree_transform = T.lprp <| mkt2 lp rp Xup
  let from_root_to_body_M = T.irootfix matmul_rev XBtoA_from_XAtoB_M (identity 6) vtree_transform
  -- let from_root_to_body_M = rootfix_work_efficient_sc matmul_rev XBtoA_from_XAtoB_M (identity 6) lp rp Xup

  let to_root_F   = map transpose  from_root_to_body_M  
  let from_root_F = map XBtoA_MtoF from_root_to_body_M 

  let fBs_root = map2 (mat_mul_vec_f64) to_root_F fBs

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 (scal_mul_vec_f64 (-1)) (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (mat_mul_vec_f64) from_root_F fJs_root 
  in map2 (vecmul) S fJs  


-- rnea with vtrees that also take an array of external forces
--  these external forces are expressed in root coordinates
-- To get an understanding of how this affects the total force see eq. 5.20 on page 95 in Roy Featherstones book
def rnea_vtree_with_f_ext [n] (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (f_ext : [n][6]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : [n]f64 =
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64)       XJ Xtree

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_M ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) 
    = (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp (zip Xup vJ)
  let vs = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
          |> map (.1)

  let S_qdd  = tabulate n (\i -> scal_mul_vec_f64 qdd[i] S[i])
  let v_cross_S_qd = tabulate n (\i -> if i == 0 then (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))
                                                 else mat_mul_vec_f64 (crm vs[i]) vJ[i])
  let as_tmp = map2 (vecadd_f64) S_qdd v_cross_S_qd 

  let vtree_as = T.lprp <| mkt2 lp rp (zip Xup as_tmp)
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
          |> map (.1) 

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let from_root_to_body_M = rootfix_work_efficient_sc matmul_rev XBtoA_from_XAtoB_M (identity 6) lp rp Xup

  let to_root_F   = map transpose  from_root_to_body_M  
  let from_root_F = map XBtoA_MtoF from_root_to_body_M 

  -- Here you add the external forces
  --  Since you transform the body forces to root coordinates the external force can just be added 
  let fBs_root = map3 (\X_to_root fbi f_xi -> 
            (X_to_root `mat_mul_vec_f64` fbi)  `vecadd_f64` (scal_mul_vec_f64 (-1) f_xi)
              ) to_root_F fBs f_ext

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 (scal_mul_vec_f64 (-1)) (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (\X_to_joint fji -> X_to_joint `mat_mul_vec_f64` fji) from_root_F fJs_root 

  in map2 (vecmul) S fJs  

def main = 
  let q = [0.8f64, 0.53f64, 0.75f64, 0.5f64, 0.91f64]
  let qd = [0.12f64, 0.85f64, 0.18f64, 0.76f64, 0.46f64]
  let qdd = [0.61f64, 0.36f64, 0.4f64, 0.73f64, 0.96f64]
  let (_, _, js, Is, Xtrees) = autoTree 5 1 0 1
  in trace <| rnea'' js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] q qd qdd [0,1,2,3,4] [9,8,7,6,5]
