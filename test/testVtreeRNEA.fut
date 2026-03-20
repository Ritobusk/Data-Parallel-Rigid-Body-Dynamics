import "../rnea"
import "../matrix_ops"
import "../spatial_ops"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"

module T = vtree
let error_tolerance = 1e-8

def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

def velocityCompare [n] (p : [n]i64) (joint_types : [n]jointT)
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
    -- let inv_cia = gauss_inv ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let vs2 = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let vs2 = map (.1) vs2

  let vs1 = loop vs' = (copy vJ) for i < (n-1) do
        let parent = p[i+1]
        in vs' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] vs'[parent]) vs'[i+1]

  in map2 (\v1 v2 -> map2 (\x y -> f64.abs (x - y) < error_tolerance) v1 v2) vs1 vs2
      |> flatten
      |> reduce (&&) true


def accelerationCompare [n] (p : [n]i64) (joint_types : [n]jointT)
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
    -- let inv_cia = gauss_inv ci.0
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

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let as2 = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let as2 = map (.1) as2

  let as1 = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
        let parent = p[i+1]
        in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) as'[i+1]



  in map2 (\a1 a2 -> map2 (\x y -> f64.abs (x - y) < error_tolerance) a1 a2) as1 as2
      |> flatten
      |> reduce (&&) true



-- ==
-- entry: velocityTest 
entry velocityTest =
  let (_, p, js, _, Is, Xtrees) = autoTree 4 2 1 1
  let lp = [0, 1, 5, 2]
  let rp = [7, 4, 6, 3]
  let t1 = velocityCompare p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0,  1] [0f64, 2, 1, 3 ] [0f64, 3, 0,  3] lp rp
  let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  let lp = [0, 1, 7, 2, 4, 8]
  let rp = [11,  6, 10, 3, 5, 9]
  let t2 = velocityCompare p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3] lp rp
  in t1 && t2

-- ==
-- entry: velocityTest 
entry accelerationTest =
  let (_, p, js, _, Is, Xtrees) = autoTree 4 2 1 1
  let lp = [0, 1, 5, 2]
  let rp = [7, 4, 6, 3]
  let t1 = accelerationCompare p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0,  1] [0f64, 2, 1, 3 ] [0f64, 3, 0,  3] lp rp
  let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  let lp = [0, 1, 7, 2, 4, 8]
  let rp = [11,  6, 10, 3, 5, 9]
  let t2 = accelerationCompare p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3] lp rp
  in t1 && t2


