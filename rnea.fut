import "matrix_ops"
import "spatial_ops"
import "treeModel"

-- Inspiration taken from https://royfeatherstone.org/spatial/v2/sourceText/ID.txt
-- As far as I understand q, qd and qdd are scalars. This might only be the case for joints with only 1-DOF
--  q, qd, qdd and tau are column vectors of length model.NB containing the joint position, velocity, acceleration and force variables, respectively.
def rnea [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd -- vector multiplication
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let (vs, as) = loop (vs', as') = (replicate n (replicate 6 0f64), replicate n (replicate 6 0f64)) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i] gravity) 
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
        let fs'' = fs' with [idx] = map2 (+) fs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fs'[idx])
        in (tau', fs'')
      else (tau', fs')
  in trace tau





def main = 
  let (_, p, js, _, Is, Xtrees) = autoTree 2 1 1 1
  in rnea p js Is Xtrees [0f64, 0, 0, 0, 0, -9.82] [0f64, 1] [0f64, 0] [0f64, 3]
