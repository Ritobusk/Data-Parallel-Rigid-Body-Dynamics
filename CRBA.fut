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
             : ([n][n]f64, [n]f64) =
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

  -- Step 2: Compute H

  let H = replicate n <| replicate n 0f64

  let I_to_root = map3 (\bXa_F aXb_M I -> bXa_F `matmul_f64`  (I `matmul_f64` aXb_M)) to_root_F from_root_to_body_M Is
  let vtree_Ics_root = T.lprp <| mkt2 lp rp I_to_root 
  let Ics_root = T.ileaffix matadd_f64 (scal_mul_mat_f64 (-1)) (replicate 6 <| replicate 6 0f64) vtree_Ics_root 
  let Ics = map3 (\bXa_F aXb_M rI -> bXa_F `matmul_f64`  (rI `matmul_f64` aXb_M)) from_root_F from_body_to_root_M Ics_root

  let depths = T.depth vtree_Ics_root
  let tree_depth = reduce i64.max 0i64 depths
  let d_sc = scan (+) 0i64  (map (+1) depths)
  let d_exsc = (rotate (-1) d_sc) with [0] = 0
  let f_arr  = 
        scatter (replicate d_sc[n-1] false) (d_exsc) (replicate n true)

  let fhs  = map2 mat_mul_vec_f64 Ics S

  let fhs' = map2 (vecmul) S fhs
  let H =  scatter_2d H (zip (iota n) (iota n)) fhs'

  let Xdown_F = map transpose Xup
  let fijs = replicate d_sc[n-1] (identity 6)

  let (_, paths) = map2 (\p d -> 
        loop (j, path) = (p, replicate tree_depth (-1i64) with [0] = 0i64) for k < d do
                let path[k] = j
                let j = parent[j]
                in (j, path)
        ) parent depths 
        |> unzip

  let paths = flatten paths |> filter (> -1i64)
  let fijs  = scatter fijs (indices paths) (map (\i -> Xdown_F[i] ) paths)
  let fijs  = scatter fijs (map (\i -> d_exsc[i]) (iota n)) 
                           (map2 (\i f -> Xdown_F[i] `matmul_f64` (diagonal f) ) (iota n) fhs)

  let fijs'  = segmented_scan matmul_rev (identity 6) f_arr fijs
             |> map get_diagonal

  -- need to scatter the (S * fijs) around H.
  -- Also I think I might have gotten confused about some index of the parent path stuff.


-- for i = 1:model.NB
--   fh = IC{i} * S{i};
--   H(i,i) = S{i}' * fh;
--   j = i;
--   while model.parent(j) > 0
--     fh = Xup{j}' * fh;
--     j = model.parent(j);
--     H(i,j) = S{j}' * fh;
--     H(j,i) = H(i,j);
--   end
-- end

  in (H, C)


