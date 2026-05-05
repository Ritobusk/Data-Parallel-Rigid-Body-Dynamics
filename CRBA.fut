import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"
import "lib/github.com/diku-dk/segmented/segmented"
import "scan_variations"

module T = vtree

def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

-- Computes the composite rigid body algorithm using vtrees.
-- The sequential implementation can be found here: https://royfeatherstone.org/spatial/v2/sourceText/HandC.txt 
-- This implementation does not do the final step of CRBA which is to solve for qdd in: H * qdd =  tau - C, for qdd
def crba_vtree [n] (parent : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) --(tau : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : ([n]f64, [n][n]f64) =
  -- Step 1: Compute the joint-space bias force: tau = ID(model, q, qd, 0)
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

  let as_tmp = tabulate n (\i -> if i == 0 then (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))
                                                 else mat_mul_vec_f64 (crm vs[i]) vJ[i])
  let vtree_as = T.lprp <| mkt2 lp rp (zip Xup as_tmp)
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
          |> map (.1) 

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let vtree_transform = T.lprp <| mkt2 lp rp Xup
  let from_root_to_body_M = T.irootfix matmul_rev XBtoA_from_XAtoB_M (identity 6) vtree_transform
  let from_body_to_root_M = map XBtoA_from_XAtoB_M from_root_to_body_M

  let to_root_F   = map transpose  from_root_to_body_M  
  let from_root_F = map XBtoA_MtoF from_root_to_body_M 

  let fBs_root = map2 (mat_mul_vec_f64) to_root_F fBs

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 (scal_mul_vec_f64 (-1)) (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (mat_mul_vec_f64) from_root_F fJs_root 
  let C = map2 (vecmul) S fJs  

  -- Step 2: Compute composite rigid bodies 

  let H = replicate n <| replicate n 0f64

  let I_to_root = map3 (\bXa_F aXb_M I -> bXa_F `matmul_f64`  (I `matmul_f64` aXb_M)) to_root_F from_root_to_body_M Is
  let vtree_Ics_root = T.lprp <| mkt2 lp rp I_to_root 
  let Ics_root = T.ileaffix matadd_f64 (scal_mul_mat_f64 (-1)) (replicate 6 <| replicate 6 0f64) vtree_Ics_root 
  let Ics = map3 (\bXa_F aXb_M rI -> bXa_F `matmul_f64`  (rI `matmul_f64` aXb_M)) from_root_F from_body_to_root_M Ics_root

  -- Step 3: Compute H

  let depths = trace <| T.depth vtree_Ics_root
  let tree_depth = reduce i64.max 0i64 depths
  let sizes = trace <| map (\x -> if x > 0 then x else 0) depths
  let d_sc = scan (+) 0i64 sizes
  let (ii1) = replicated_iota sizes 
  let ii1 = sized d_sc[n-1] ii1
  let (ii1) = trace (ii1)

  let fhs  = map2 mat_mul_vec_f64 Ics S

  let fhs' = map2 (vecmul) S fhs
  let H =  scatter_2d H (zip (iota n) (iota n)) fhs'

  let parent = trace parent
  let (_, paths) = trace <| map2 (\p d -> 
        loop (j, path) = (p, replicate tree_depth (-1i64) ) for k < (d) do
                let path = path with [k] = j 
                let j = parent[j]
                in (j, path)
        ) parent depths 
        |> unzip
  let paths = flatten paths |> filter (> -1i64) |> sized (d_sc[n-1])

  let tmp = map2 (\i j ->
          let f = from_root_F[j] `mat_mul_vec_f64` (to_root_F[i] `mat_mul_vec_f64` fhs[i])
          let fh = S[j] `vecmul` f
          in fh
      ) ii1 paths

  let H =  scatter_2d H (zip ii1 paths) tmp
  let H =  scatter_2d H (zip paths ii1) tmp

  in (C, H)

def crba_vtree' [n] [nd] [ndd] (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) --(tau : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             (paths : [nd]i64) (p_ii1 : [ndd]i64)
             : ([n]f64, [n][n]f64) =
  -- Step 1: Compute the joint-space bias force: tau = ID(model, q, qd, 0)
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

  let as_tmp = tabulate n (\i -> if i == 0 then (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))
                                                 else mat_mul_vec_f64 (crm vs[i]) vJ[i])
  let vtree_as = T.lprp <| mkt2 lp rp (zip Xup as_tmp)
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_as
          |> map (.1) 

  let fBs = tabulate n 
          (\i -> Is[i] `mat_mul_vec_f64` as[i] `vecadd_f64` ((crf vs[i]) `matmul_f64 ` Is[i] `mat_mul_vec_f64` vs[i]))

  let vtree_transform = T.lprp <| mkt2 lp rp Xup
  let from_root_to_body_M = T.irootfix matmul_rev XBtoA_from_XAtoB_M (identity 6) vtree_transform
  let from_body_to_root_M = map XBtoA_from_XAtoB_M from_root_to_body_M

  let to_root_F   = map transpose  from_root_to_body_M  
  let from_root_F = map XBtoA_MtoF from_root_to_body_M 

  let fBs_root = map2 (mat_mul_vec_f64) to_root_F fBs

  let vtree_fjs_root = T.lprp <| mkt2 lp rp fBs_root
  let fJs_root = T.ileaffix vecadd_f64 (scal_mul_vec_f64 (-1)) (replicate 6 0f64) vtree_fjs_root

  let fJs = map2 (mat_mul_vec_f64) from_root_F fJs_root 
  let C = map2 (vecmul) S fJs  

  -- Step 2: Compute composite rigid bodies 

  let H = replicate n <| replicate n 0f64

  let I_to_root = map3 (\bXa_F aXb_M I -> bXa_F `matmul_f64`  (I `matmul_f64` aXb_M)) to_root_F from_root_to_body_M Is
  let vtree_Ics_root = T.lprp <| mkt2 lp rp I_to_root 
  let Ics_root = T.ileaffix matadd_f64 (scal_mul_mat_f64 (-1)) (replicate 6 <| replicate 6 0f64) vtree_Ics_root 
  let Ics = map3 (\bXa_F aXb_M rI -> bXa_F `matmul_f64`  (rI `matmul_f64` aXb_M)) from_root_F from_body_to_root_M Ics_root

  -- Step 3: Compute H
  let fhs  = map2 mat_mul_vec_f64 Ics S
  let Hii = map2 (vecmul) S fhs
  let H =  scatter_2d H (zip (iota n) (iota n)) Hii
  let p_ii1 = sized nd p_ii1

  let Hij = map2 (\i j ->
          let f = from_root_F[j] `mat_mul_vec_f64` (to_root_F[i] `mat_mul_vec_f64` fhs[i])
          in S[j] `vecmul` f
      ) p_ii1 paths

  let H =  scatter_2d H (zip p_ii1 paths) Hij
  let H =  scatter_2d H (zip paths p_ii1) Hij
  in (C, H)



def main = 
  let q = [0.8f64, 0.53f64, 0.75f64, 0.5f64, 0.91f64]
  let qd = [0.12f64, 0.85f64, 0.18f64, 0.76f64, 0.46f64]
  let (_, _, js, Is, Xtrees,lp,rp, paths, p_ii1) = autoVTree 5 1 0 1
  in trace <| crba_vtree' js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] q qd lp rp paths p_ii1
