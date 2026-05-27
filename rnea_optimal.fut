import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "vtree_with_work_efficient_scan"

module T = vtree

def mkt 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

def rnea_vtree_optimized_ds [n] (joint_types : [n]jointT)
             (Is : [n]I_Compact) (Xtree : [n]X_Compact) 
             (gravity : mv)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : [n]f64 =

  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) joint_types q
  let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtree

  let vtree_transform = T.lprp <| mkt lp rp Xup
  let transformation_tree = T.irootfix_b transform_XX_rev X_inv (copy transform_identity) vtree_transform 64i64

  let (vJ, aJ) = map4 (\qdi qddi si Xroot -> (Xroot `Xm_inv` (qdi `scal_mul_mv` si),
                                              Xroot `Xm_inv` (qddi `scal_mul_mv` si))
                       ) qd qdd S transformation_tree
                |> unzip
  
  let vs = irootfix_vector_add lp rp (map mv_to_6d vJ)
            |> map d6_to_mv 

  let as_tmp = tabulate n (\i -> if i == 0 then aJ[i] `mv_add` (scal_mul_mv (-1) gravity)
                                           else aJ[i] `mv_add` (vs[i] `mv_cross_mv` vJ[i]))

  let as_root = irootfix_vector_add  lp rp (map mv_to_6d as_tmp)
                |> map d6_to_mv

  let vs = map2 Xm transformation_tree vs
  let as = map2 Xm transformation_tree as_root

  let fBs_root = tabulate n 
          (\i -> transformation_tree[i] `Xf_inv` 
                     (Is[i] `IC_mul_mv` as[i] `fv_add` (vs[i] `mv_cross_fv` (Is[i] `IC_mul_mv` vs[i])))
          )
      |> map fv_to_6d

  let fJs_root = ileaffix_vector_add  lp rp fBs_root
      |> map d6_to_fv 

  in map3 (\frt fji si -> si `scalar_prod` (frt `Xf` fji))  transformation_tree fJs_root S

def inv_op_C (ci : (X_Compact, mv)) : (X_Compact, mv) =
      let inv_cia = transform_inv ci.0
      let inv_cib = scal_mul_mv (-1) (inv_cia `Xm` ci.1)
      in (inv_cia, inv_cib)
def operator_C (si : (X_Compact, mv)) (ci : (X_Compact, mv)) : (X_Compact, mv) =
      (ci.0 `transform_XX` si.0,    (ci.0 `Xm` si.1) `mv_add` ci.1)

-- Using the bullet op for transformation tree and velocities
def rnea_vtree_optimized_ds' [n] (joint_types : [n]jointT)
             (Is : [n]I_Compact) (Xtree : [n]X_Compact) 
             (gravity : mv)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             : [n]f64 =

  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) joint_types q
  let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtree
  let (vJ, aJ) = map3 (\qdi qddi si -> (qdi `scal_mul_mv` si, qddi `scal_mul_mv` si)) qd qdd S
                |> unzip
  let vtree_transform = T.lprp <| mkt lp rp (zip Xup vJ)
  let (transformation_tree, vs) = T.irootfix_b  operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) vtree_transform 64i64
    |> unzip

  let as_tmp = tabulate n (\i -> if i == 0 then aJ[i] `mv_add` ( Xup[0] `Xm` ((-1) `scal_mul_mv` gravity))
                                           else aJ[i] `mv_add` (vs[i] `mv_cross_mv` vJ[i])
                          )
            |> map2 Xm_inv transformation_tree 

  let as_root = irootfix_vector_add  lp rp (map mv_to_6d as_tmp)
                |> map d6_to_mv

  let as = map2 Xm transformation_tree as_root

  let fBs_root = tabulate n 
          (\i -> transformation_tree[i] `Xf_inv` 
                     (Is[i] `IC_mul_mv` as[i] `fv_add` (vs[i] `mv_cross_fv` (Is[i] `IC_mul_mv` vs[i])))
          )
      |> map fv_to_6d

  let fJs_root = ileaffix_vector_add  lp rp fBs_root
      |> map d6_to_fv 

  in map3 (\frt fji si -> si `scalar_prod` (frt `Xf` fji))  transformation_tree fJs_root S




def rnea_optimized_ds_seq [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n]I_Compact) (Xtree : [n]X_Compact) 
             (gravity : mv) (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             : [n]f64 =
  let (XJ, S) = unzip <| map2 (jcalcC) joint_types q 
  let vJ      = map2 (scal_mul_mv) qd  S 
  let Xup     = map2 (transform_XX) XJ Xtree

  let empty_mv = {w = [0,0,0f64], v_O = [0,0,0f64]}
  let (vs, as) = loop (vs', as') = (replicate n empty_mv, replicate n empty_mv) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = (mv_add) 
                                   (Xup[i] `Xm`  ( (-1) `scal_mul_mv` gravity  )) 
                                   (qdd[i] `scal_mul_mv` S[i])
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = mv_add (Xup[i] `Xm`  (copy vs'[parent]))  vJ[i]

        let as' = as' with [i] = mv_add ( Xup[i] `Xm` (copy as'[parent]))   
                                           ( qdd[i] `scal_mul_mv` S[i])
                                |> mv_add (  vs'[i] `mv_cross_mv` vJ[i])
        in (vs', as')

  let fBs = tabulate n 
          (\i -> Is[i] `IC_mul_mv` as[i] `fv_add` ( vs[i] `mv_cross_fv` (Is[i] `IC_mul_mv` vs[i])))

  let (tau, _) = loop (tau', fs') = (replicate n 0f64, fBs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = scalar_prod S[idx]   fs'[idx] 
    in 
      if idx > 0 then
        let tmp =  fv_add (fs'[parent]) (Xup[idx] `Xf_inv` (fs'[idx]))
        let fs'' = fs' with [parent] = copy tmp
        in (tau', fs'')
      else (tau', fs')
  in tau


