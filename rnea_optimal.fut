import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "vtree_with_work_efficient_scan"

module T = vtree

def rootfix_vector_add [n]   (lp : [n]i64) (rp : [n]i64) (data : [n][6]f64) : [n][6]f64 =
    let I = replicate (2 * n) (replicate 6 0f64)
    let L = scatter I lp data
    let R = scatter L rp (map ((scal_mul_vec_f64 (-1))) data)
    let S = transpose R
            |> map (exscan (+) 0f64)
            |> transpose
    in map (\i -> S[i]) lp

def ileaffix_vector_add 'a [n] (lp : [n]i64) (rp : [n]i64) (data : [n][6]f64) : [n][6]f64 =
    let I = replicate (2 * n) (replicate 6 0f64)
    let L = scatter I lp data
    let S = transpose L
            |> map (exscan (+) 0f64)
            |> transpose
    let Rv = map (\i -> S[i]) rp
    let Lv = map (\i -> scal_mul_vec_f64 (-1) S[i]) lp
    in map2 vecadd_f64 Rv Lv

def irootfix_vector_add 'a [n] (lp: [n]i64) (rp: [n]i64) (data : [n][6]f64) : [n][6]f64 =
    map2 vecadd_f64 (rootfix_vector_add lp rp data) data

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

