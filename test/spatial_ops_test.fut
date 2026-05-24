import "../treeModel"
import "../spatial_ops"
import "../matrix_ops"
import "../lib/github.com/diku-dk/cpprandom/random"
import "../vtree_with_work_efficient_scan"
import "../rnea_optimal"

module T = vtree

module mktest (dist: rng_distribution) = {
  module engine = dist.engine
  module num = dist.num

  def test (x: i32) (n: i64) (d: dist.distribution) =
    let rng = engine.rng_from_seed [x]
    let (rng, _) = dist.rand d rng
    let rngs = engine.split_rng n rng
    let (rngs', qs) = unzip (map (dist.rand d) rngs)
    let (rngs'', qds) = unzip (map (dist.rand d) rngs')
    let (_, qdds) = unzip (map (dist.rand d) rngs'')
    in (qs, qds, qdds) 
}

module test_f32_rand_m =
  mktest (uniform_real_distribution f32 minstd_rand)

let error_tolerance = 1e-7

def compare_matrices [n][l][m] (X1: [n][l][m]f64) (X2: [n][l][m]f64) : bool =
  map2 (\x1 x2 -> map2 (\r r' -> map2 (\x y -> f64.abs (x - y) < error_tolerance) r r') x1 x2) X1 X2
      |> flatten |> flatten
      |> reduce (&&) true

def compare_vectors [n][m] (V1: [n][m]f64) (V2: [n][m]f64) : bool =
  map2 (\v1 v2 ->  map2 (\x y -> f64.abs (x - y) < error_tolerance) v1 v2) V1 V2
      |> flatten 
      |> reduce (&&) true

def mkt 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

def test_autoTree  (n : i64) :   
    ([n]jointT, [n][6][6]f64, [n][6][6]f64, [n]i64, [n]i64, [n]f64, [n]f64, [n]f64 ) =
    let (q, qd, qdd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    let qdd = map f64.f32 qdd
    let (_, p, js, Is, Xtrees) = autoTree n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    in (js, Is, Xtrees, lp, rp, q, qd, qdd)

def test_autoTreeC  (n : i64) :   
    ([n]jointT, [n]I_Compact, [n]X_Compact, [n]i64, [n]i64, [n]f64, [n]f64, [n]f64 ) =
    let (q, qd, qdd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    let qdd = map f64.f32 qdd
    let (_, p, js, Is, Xtrees) = autoTreeC n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    in (js, Is, Xtrees, lp, rp, q, qd, qdd)

-- Tests the optimized/compact data structures against the normal 6d vectors and 6x6 matrices of spatial vector algebra

def I_multiplication (n : i64) : bool =
  let (js, Is, _,_,_, q, _, _) = test_autoTree n

  let (_, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js q 
  let fS     = map2 (\is s -> mat_mul_vec_f64 is s) Is S

  let (jsC, IsC, _,_,_, _, _, _) = test_autoTreeC n
  let (_, SC) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) jsC q 
  let fS'     = map2 (\is s -> IC_mul_mv is s) IsC SC
  let fS' = map (fv_to_6d) fS'

  in compare_vectors fS fS' 


def transform_multiplication (n : i64) : bool =
  let (js, _, Xtrees,_,_, q, _, _) = test_autoTree n

  let (XJ, _) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js q 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtrees

  let (jsC, _, XtreesC,_,_,_, _, _) = test_autoTreeC n
  let (XJC, _) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) jsC q 
  let XupC     = map2 (\xj xtree -> transform_XX xj xtree) XJC XtreesC
  let XupC' = map (from_XC_to_aXb) XupC

  in compare_matrices Xup XupC' 

def rootfix_vj (n : i64) : [3]bool =
  let (joint_types, Is, Xtree, lp, rp, q, qd, qdd) = test_autoTree n
  let gravity = [0,0,0,0,0,-9.81f64]
  let gravity_mv = d6_to_mv gravity

  let (Xup, S) = unzip <| map3 (\Xti jti qi -> 
                    let (Xji, si) = jcalc jti qi
                    in (Xji `matmul_f64` Xti, si)) Xtree joint_types q 

  let vtree_transform = T.lprp <| mkt lp rp Xup
  let transformation_tree = T.irootfix matmul_rev XBtoA_from_XAtoB_M (identity 6) vtree_transform
    |> trace

  let to_root_F   = map transpose  transformation_tree  
  let from_root_F = map XBtoA_MtoF transformation_tree 
  let to_root_M =   map (XBtoA_FtoM) to_root_F 

  let (vJ, aJ) = map4 (\qdi qddi si Xroot -> (Xroot `mat_mul_vec_f64` (qdi `scal_mul_vec_f64` si),
                                              Xroot `mat_mul_vec_f64` (qddi `scal_mul_vec_f64` si))
                       ) qd qdd S to_root_M
                |> unzip
  let vs = irootfix_vector_add  lp rp vJ

  let as_tmp = tabulate n (\i -> if i == 0 then aJ[i] `vecadd_f64` ( (map (\x -> -1 * x) gravity))
                                           else aJ[i] `vecadd_f64` (mat_mul_vec_f64 (crm vs[i]) vJ[i]))

  let as_root = irootfix_vector_add  lp rp as_tmp

  let (joint_typesC, IsC, XtreeC, _,_,_,_,_) = test_autoTreeC n
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) joint_typesC q
  let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ XtreeC

  let vtree_transformC = T.lprp <| mkt lp rp Xup
  let transformation_treeC = T.irootfix_b transform_XX_rev X_inv (copy transform_identity) vtree_transformC 64i64

  let res1 = compare_matrices transformation_tree (map from_XC_to_aXb transformation_treeC)

  let (vJ, aJ) = map4 (\qdi qddi si Xroot -> (Xroot `Xm_inv` (qdi `scal_mul_mv` si),
                                              Xroot `Xm_inv` (qddi `scal_mul_mv` si))
                       ) qd qdd S transformation_treeC
                |> unzip
  
  let vsC = irootfix_vector_add lp rp (map mv_to_6d vJ)
  let res2 = compare_vectors vs vsC
  let vsC = map d6_to_mv  vsC
  

  let as_tmp = tabulate n (\i -> if i == 0 then aJ[i] `mv_add` (scal_mul_mv (-1) gravity_mv)
                                           else aJ[i] `mv_add` (vsC[i] `mv_cross_mv` vJ[i]))

  let as_rootC = irootfix_vector_add  lp rp (map mv_to_6d as_tmp)
  let res3 = compare_vectors as_root as_rootC
  let as_rootC = map d6_to_mv as_rootC 


  in [res1, res2, res3]



-- ==
-- entry: plucker_transform 
-- input {4i64} output {true}
-- input {8i64} output {true}
entry plucker_transform (n : i64) =
  transform_multiplication n

-- ==
-- entry: IC_mul 
-- input {4i64} output {true}
-- input {8i64} output {true}
entry IC_mul (n : i64) =
  I_multiplication n

-- ==
-- entry: mat_mat_rootfix 
-- input {4i64} output {[true, true, true]}
-- input {8i64} output {[true, true, true]}
entry mat_mat_rootfix (n : i64) =
  rootfix_vj n

def main =
  rootfix_vj 4
