import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "vtree_with_work_efficient_scan"
import "lib/github.com/diku-dk/segmented/segmented"
import "CRBA"

module T = vtree

def exscan f ne xs =
    map2 (\i x -> if i == 0 then ne else x)
     (indices xs)
     (rotate (-1) (scan f ne xs))

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

def crba_seq_optimized_ds [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n]I_Compact) (Xtree : [n]X_Compact) 
             (gravity : mv)
             (q : [n]f64) (qd : [n]f64)
            : ([n]f64, [n][n]f64) =

  let (XJ, S) = unzip <| map2 (jcalcC) joint_types q 
  let vJ      = map2 (scal_mul_mv) qd  S 
  let Xup     = map2 (transform_XX) XJ Xtree

  let empty_mv = {w = [0,0,0f64], v_O = [0,0,0f64]}
  let (vs, as) = loop (vs', as') = (replicate n empty_mv, replicate n empty_mv) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = Xm Xup[i] ( (-1) `scal_mul_mv` gravity  )
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = mv_add (Xup[i] `Xm`  (copy vs'[parent]))  vJ[i]

        let as' = as' with [i] = mv_add ( Xup[i] `Xm` (copy as'[parent]))   
                                 (vs'[i] `mv_cross_mv` vJ[i])
        in (vs', as')

  let fBs = tabulate n 
          (\i -> Is[i] `IC_mul_mv` as[i] `fv_add` ( vs[i] `mv_cross_fv` (Is[i] `IC_mul_mv` vs[i])))

  let (C, _) = loop (tau', fs') = (replicate n 0f64, fBs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = scalar_prod S[idx]   fs'[idx] 
    in 
      if idx > 0 then
        let fs'' = fs' with [parent] = (fv_add) (copy fs'[parent]) (Xup[idx] `Xf_inv` (copy fs'[idx]))
        in (tau', fs'')
      else (tau', fs')

  -- Step 2: Compute composite rigid bodies 
  let Ics = loop IC = (copy Is) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    in 
        if idx > 0 then
          let tmp = (copy IC[parent]) `IC_add` (transform_IC_inv Xup[idx] IC[idx]) 
          let ic = IC with [parent] = copy tmp
          in ic
        else IC

  let p' = map (\i -> if i == 0 then -1 else p[i]) (iota n)

  let H'' = loop H = (replicate n <| replicate n 0f64) for i < n do
    let fh = Ics[i] `IC_mul_mv` S[i]
    let H_tmp  = H with [i,i] = scalar_prod S[i] fh
    let (H'',_,_) =  loop (h, j, f) = (H_tmp, (copy i), (copy fh)) while p'[j] >= 0 do
        let fh' = (Xup[j]) `Xf_inv` f
        let j = p'[j]
        let tmp = scalar_prod S[j] fh'
        let h' = h with [i,j] = tmp
        let h'' = h' with [j,i] = (copy h'[i,j])
        in (h'', j, fh')
    in H''
  in (C, H'')


 

def crba_vtree_optimized_ds [n] [nd] [ndd] (joint_types : [n]jointT)
             (Is : [n]I_Compact) (Xtree : [n]X_Compact) 
             (gravity : mv)
             (q : [n]f64) (qd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) 
             (paths : [nd]i64) (p_ii1 : [ndd]i64)
             : ([n]f64, [n][n]f64) =
  -- Step 1: Compute the joint-space bias force: tau = ID(model, q, qd, 0)
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) joint_types q
  let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtree

  let vtree_transform = T.lprp <| mkt lp rp Xup
  let transformation_tree = T.irootfix_b transform_XX_rev X_inv (copy transform_identity) vtree_transform 64i64

  let vJ = map3 (\qdi si Xroot -> (Xroot `Xm_inv` (qdi `scal_mul_mv` si)) ) qd S transformation_tree
  
  let vs = irootfix_vector_add lp rp (map mv_to_6d vJ)
            |> map d6_to_mv 

  let as_tmp = tabulate n (\i -> if i == 0 then (scal_mul_mv (-1) gravity)
                                           else (vs[i] `mv_cross_mv` vJ[i]))

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

  let C =  map3 (\frt fji si -> si `scalar_prod` (frt `Xf` fji))  transformation_tree fJs_root S

  -- Step 2: Compute composite rigid bodies 
  let I_to_root = map2 (transform_IC_inv) transformation_tree Is
  let vtree_Ics_root = T.lprp <| mkt lp rp I_to_root 
  let Ics_root = T.ileaffix IC_add IC_inv (copy IC_identity) vtree_Ics_root 
  let Ics = map2 (transform_IC) transformation_tree Ics_root

  -- Step 3: Compute H
  let H = replicate n <| replicate n 0f64

  let fhs  = map2 IC_mul_mv Ics S
  let Hii = map2 (scalar_prod) S fhs
  let H =  scatter_2d H (zip (iota n) (iota n)) Hii
  let p_ii1 = sized nd p_ii1

  let Hij = map2 (\i j ->
          let f = transformation_tree[j] `Xf` (transformation_tree[i] `Xf_inv` fhs[i])
          in S[j] `scalar_prod` f
      ) p_ii1 paths

  let H =  scatter_2d H (zip p_ii1 paths) Hij
  let H =  scatter_2d H (zip paths p_ii1) Hij
  in (C, H)

def main = 
  let q = [0.8f64, 0.53f64, 0.75f64, 0.5f64, 0.91f64]
  let qd = [0.12f64, 0.85f64, 0.18f64, 0.76f64, 0.46f64]
  let (_, p, js, Is, Xtrees,lp,rp, paths, p_ii1) = autoVTree 5 1 0 1
  let (_, _, _, IsC, XtreesC,_,_,_,_) = autoVTreeC 5 1 0 1
  let gravity = {w = [0,0,0f64], v_O = [0,0, -9.81f64]}
  let (_,_) = crba_seq_optimized_ds p js IsC XtreesC gravity q qd 
  in  crba_seq p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] q qd 
