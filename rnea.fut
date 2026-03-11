import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"

module T = vtree

def mkt 'a [n] (ps:[n]i64) (ds:[n]a) : [n]{parent:i64,data:a} =
    map2 (\p d -> {parent=p,data=d}) ps ds

-- Inspiration taken from https://royfeatherstone.org/spatial/v2/sourceText/ID.txt
-- As far as I understand q, qd and qdd are scalars. This might only be the case for joints with only 1-DOF
--  q, qd, qdd and tau are column vectors of length model.NB containing the joint position, velocity, acceleration and force variables, respectively.
def rnea [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let (vs, as) = loop (vs', as') = (replicate n (replicate 6 0f64), replicate n (replicate 6 0f64)) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i] (map (\x -> -1 * x) gravity)) 
                                          (map (\s -> s * qdd[i]) S[i])
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i] vs'[parent]) 
                                           vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i] as'[parent]) -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
                                          (map (\s -> s * qdd[i]) S[i])
                                |> map2 (+) (mat_mul_vec_f64 (crm vs'[i]) vJ[i])
        in (vs', as')

  let fs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let (tau, _) = loop (tau', fs') = (replicate n 0f64, fs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = vec_mul_vec S[idx] fs'[idx] 
    in 
      if idx > 0 then
        let fs'' = fs' with [parent] = map2 (+) fs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fs'[idx])
        in (tau', fs'')
      else (tau', fs')
  in trace tau


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let vs = loop vs' = (copy vJ) for i < (n-1) do
        let parent = p[i+1]
        in vs' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] vs'[parent]) vs'[i+1]

  let as_tmp = map2 (\S_qdd v_cross_S_qd -> map2 (+) S_qdd v_cross_S_qd)
                        (map (\i -> map (\s -> s * qdd[i]) S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then replicate 6 0f64
                          else mat_mul_vec_f64 (crm vs[i]) vJ[i]) (iota n))
  let as_tmp = as_tmp with [0] = map2 (+) as_tmp[0] (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))

  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
        let parent = p[i+1]
        in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) as'[i+1]

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let fJs = loop fJs' = (fBs) for i < n -1 do
    let idx = n - (i+1)
    let parent = p[idx]
    in fJs' with [parent] = map2 (+) fJs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fJs'[idx])
  in trace <| map2 (\s f -> vec_mul_vec s f) S fJs  

-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea'' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = trace <| map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = trace <| map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let test = trace <| map (\x -> (gauss_jordan x) ) Xup
  
  let vtree_vs = T.mk_preorder <| mkt p <| zip vJ <| iota n 
  let operator (p_vec_id : ([6]f64, i64)) (vec_id : ([6]f64, i64)) = (map2 (+) (mat_mul_vec_f64 Xup[p_vec_id.1] p_vec_id.0) vec_id.0, vec_id.1) 
  let inv_operator (vec_id : ([6]f64, i64)) =
      let inv_mat = gauss_jordan Xup[vec_id.1]
      in (mat_mul_vec_f64 inv_mat vec_id.0, vec_id.1) 

  let vs2 = T.ileaffix operator inv_operator (replicate 6 0f64, 0) vtree_vs

  let vs = trace <| loop vs' = (copy vJ) for i < (n-1) do
        let parent = p[i+1]
        in vs' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] vs'[parent]) vs'[i+1]

  let as_tmp = map2 (\S_qdd v_cross_S_qd -> map2 (+) S_qdd v_cross_S_qd)
                        (map (\i -> map (\s -> s * qdd[i]) S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then replicate 6 0f64
                          else mat_mul_vec_f64 (crm vs[i]) vJ[i]) (iota n))
  let as_tmp = as_tmp with [0] = map2 (+) as_tmp[0] (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))

  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
        let parent = p[i+1]
        in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) as'[i+1]

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let fJs = loop fJs' = (fBs) for i < n -1 do
    let idx = n - (i+1)
    let parent = p[idx]
    in fJs' with [parent] = map2 (+) fJs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fJs'[idx])
  in trace test
  --in trace <| map2 (\s f -> vec_mul_vec s f) S fJs  


def main = 
  let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3]
