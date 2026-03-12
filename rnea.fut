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
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (vec_mul_scalar_f64) S qd 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let (vs, as) = loop (vs', as') = (replicate n (replicate 6 0f64), replicate n (replicate 6 0f64)) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i]  ( gravity `vec_mul_scalar_f64` (-1) )) 
                                   (S[i] `vec_mul_scalar_f64` qdd[i])
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = map2 (+) (Xup[i] `mat_mul_vec_f64`  vs'[parent])  vJ[i]

        let as' = as' with [i] = map2 (+) ( Xup[i] `mat_mul_vec_f64` as'[parent])   
                                           (S[i] `vec_mul_scalar_f64` qdd[i])
                                |> vecadd_f64 ( (crm vs'[i]) `mat_mul_vec_f64` vJ[i])
        in (vs', as')

  let fBs = map (\i -> 
              (mat_mul_vec_f64 Is[i] as[i]) `vecadd_f64` (((crf vs[i]) `matmul_f64` Is[i]) `mat_mul_vec_f64` vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let (tau, _) = loop (tau', fs') = (replicate n 0f64, fBs) for i < n do
    let idx = n - (i+1)
    let parent = p[idx]
    let tau' = tau' with [idx] = vecmul S[idx]  fs'[idx] 
    in 
      if idx > 0 then
        let fs'' = fs' with [parent] = map2 (+) fs'[parent] ((transpose Xup[idx]) `mat_mul_vec_f64` fs'[idx])
        in (tau', fs'')
      else (tau', fs')
  in trace tau


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (vec_mul_scalar_f64) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let vs = loop vs' = (copy vJ) for i < (n-1) do
        let parent = p[i+1]
        in vs' with [i+1] = map2 (+) (Xup[i+1] `mat_mul_vec_f64` vs'[parent]) vs'[i+1]

  let as_tmp = map2 (vecadd_f64)
                        (map (\i ->  S[i] `vec_mul_scalar_f64` qdd[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity)
                          else  (crm vs[i]) `mat_mul_vec_f64` vJ[i]) (iota n))

  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
        let parent = p[i+1]
        in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) as'[i+1]

  let fBs = map (\i -> 
              (Is[i] `mat_mul_vec_f64` as[i]) `vecadd_f64` (((crf vs[i]) `matmul_f64` Is[i]) `mat_mul_vec_f64` vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let fJs = loop fJs' = (fBs) for i < n -1 do
    let idx = n - (i+1)
    let parent = p[idx]
    in fJs' with [parent] = map2 (+) fJs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fJs'[idx])
  in trace <| map2 (\s f -> s `vecmul`  f) S fJs  


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea'' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = trace <| map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree
  -- let Xup     = Xup ++ ([tabulate_2d 6 6 (\r c -> if r == c then 1f64 else 0f64)])
  
  let combine_mat_vec mat v = tabulate_2d 7 7 
                    (\r c -> if c < 6 && r < 6 then mat[r][c]
                             else if r < 6 && c == 6 then v[r]  
                             else 0f64
                    )


  let A = map2 (\x vj -> combine_mat_vec x vj) Xup vJ


  -- let vtree_vs = T.mk_preorder <| mkt p <| zip vJ <| iota n 
  -- let operator (p_vec_id : ([6]f64, i64)) (vec_id : ([6]f64, i64)) = (map2 (+) (mat_mul_vec_f64 Xup[vec_id.1] p_vec_id.0) vec_id.0, vec_id.1) 
  -- let inv_operator (vec_id : ([6]f64, i64)) =
  --     let inv_mat = gauss_jordan Xup[vec_id.1]
  --     in (mat_mul_vec_f64 inv_mat vec_id.0, vec_id.1) 

  let vtree_vs = T.mk_preorder <| mkt p A
  let operator (p_a : [7][7]f64) (a : [7][7]f64) = 
      let new_vec = mat_mul_vec_f64 a[:6, :6] p_a[:6, 6]
      let new_vec = map2 (+) new_vec a[:6, 6]
      in combine_mat_vec a[:6, :6] new_vec
  let inv_operator (a : [7][7]f64) =
      let inv_mat = gauss_jordan a[:6, :6]
      let new_vec = map (\x -> x * (-1)) a[:6, 6]
      in combine_mat_vec a[:6, :6] new_vec
  -- let inv_operator (a : [7][7]f64) =
  --     let inv_mat = gauss_jordan a[:6, :6]
  --     let new_vec = map2 (-) mat_mul_vec_f64 inv_mat a[:6, 6]
  --     in combine_mat_vec a[:6, :6] new_vec

  -- let vs2 = T.irootfix operator inv_operator (replicate 6 0f64, 0) vtree_vs
  let vs2 = T.irootfix operator inv_operator (identity 7) vtree_vs

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
  in trace vs2
  --in trace <| map2 (\s f -> vecmul s f) S fJs  


def main = 
  let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  in rnea' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3]


-- Ide: 1. Lav en data struktur, så man ved hvilke links er i hvilke dybder af det kinematiske træ
--      2. Kør sequentielt gennem dybderne
--      3. Ved hver dybde udregn velocity eller acceleration 

-- Tanke: Tror ikke vs og as kan beregnes med vtrees, da deres operatorer ikke synes at være assosiative.
