import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"

module T = vtree

def mkt 'a [n] (ps:[n]i64) (ds:[n]a) : [n]{parent:i64,data:a} =
    map2 (\p d -> {parent=p,data=d}) ps ds
def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

-- Inspiration taken from https://royfeatherstone.org/spatial/v2/sourceText/ID.txt
-- As far as I understand q, qd and qdd are scalars. This might only be the case for joints with only 1-DOF
--  q, qd, qdd and tau are column vectors of length model.NB containing the joint position, velocity, acceleration and force variables, respectively.
def rnea [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64) = 
  let (XJ, S) = unzip <| map2 (jcalc) joint_types q 
  let vJ      = map2 (scal_mul_vec_f64) qd  S 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let (vs, as) = loop (vs', as') = (replicate n (replicate 6 0f64), replicate n (replicate 6 0f64)) for i < n do
    if i == 0 then 
        let vs' = vs' with [i] = vJ[i]
        let as' = as' with [i] = map2 (+) (mat_mul_vec_f64 Xup[i]  ( (-1) `scal_mul_vec_f64` gravity  )) 
                                   (qdd[i] `scal_mul_vec_f64` S[i])
        in (vs', as')
    else 
        let parent = p[i]
        let vs' = vs' with [i] = map2 (+) (Xup[i] `mat_mul_vec_f64`  vs'[parent])  vJ[i]

        let as' = as' with [i] = map2 (+) ( Xup[i] `mat_mul_vec_f64` as'[parent])   
                                           ( qdd[i]`scal_mul_vec_f64` S[i])
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
  let vJ      = map2 (scal_mul_vec_f64) qd S 
  let Xup     = map2 (matmul_f64) XJ Xtree

  let vs = loop vs' = (copy vJ) for i < (n-1) do
      let parent = p[i+1]
      in vs' with [i+1] = map2 (+) (Xup[i+1] `mat_mul_vec_f64` vs'[parent]) 
                                    vs'[i+1]

  let as_tmp = map2 (vecadd_f64)
                        (map (\i -> qdd[i]  `scal_mul_vec_f64` S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity)
                          else  (crm vs[i]) `mat_mul_vec_f64` vJ[i]) (iota n))

  let as = loop as' = as_tmp for i < n-1 do  -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
      let parent = p[i+1]
      in as' with [i+1] = map2 (+) (mat_mul_vec_f64 Xup[i+1] as'[parent]) 
                                     as'[i+1]

  let fBs = map (\i -> 
              (Is[i] `mat_mul_vec_f64` as[i]) `vecadd_f64` (((crf vs[i]) `matmul_f64` Is[i]) `mat_mul_vec_f64` vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let fJs = loop fJs' = (fBs) for i < n -1 do
      let idx = n - (i+1)
      let parent = p[idx]
      let fJs'' = fJs' with [parent] = map2 (+) (copy fJs'[parent])
                                    <| mat_mul_vec_f64 (transpose Xup[idx]) (copy fJs'[idx])
      in fJs''

  in trace <| map2 (\s f -> s `vecmul`  f) S fJs  


-- Same as above implementation but it uses a bit more parallelism such that only
--  the computations that v-trees are necessary for are exposed.
def rnea'' [n] (p : [n]i64) (joint_types : [n]jointT)
             (Is : [n][6][6]f64) (Xtree : [n][6][6]f64) 
             (gravity : [6]f64)
             (q : [n]f64) (qd : [n]f64) (qdd : [n]f64)
             (lp : [n]i64) (rp : [n]i64) =
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) joint_types q 
  let vJ      = map2 (\s v -> map (\x -> x * v) s) S qd 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtree

  let Cs = zip Xup vJ

  let p = trace p

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_M ci.0
    -- let inv_cia = gauss_inv ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)
  
  let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let vs = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  -- let vtree_vs = T.mk_preorder <| mkt p Cs
  -- let vs2 = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let vs = map (.1) vs

  -- as_tmp = S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 
  let as_tmp = map2 (\S_qdd v_cross_S_qd -> map2 (+) S_qdd v_cross_S_qd)
                        (map (\i -> map (\s -> s * qdd[i]) S[i]) (iota n))
                        (map (\i -> 
                          if   i == 0 then replicate 6 0f64
                          else mat_mul_vec_f64 (crm vs[i]) vJ[i]) (iota n))
  let as_tmp = as_tmp with [0] = map2 (+) as_tmp[0] (mat_mul_vec_f64 Xup[0] (map (\x -> -1 * x) gravity))

  let Cs = zip Xup as_tmp

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let as = T.irootfix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let as = map (.1) as -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 

  let fBs = map (\i -> 
              map2 (+) (mat_mul_vec_f64 Is[i] as[i]) (mat_mul_vec_f64 (matmul_f64 (crf vs[i]) Is[i]) vs[i])
              ) (iota n) 

  -- Here you should add the external forces. This is not yet implemented!
  --   f = apply_external_forces( model.parent, Xup, f, f_ext );

  let Cs = zip (map (transpose) Xup) fBs

  let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    let inv_cia = XBtoA_from_XAtoB_F ci.0
    let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    in (inv_cia, inv_cib)

  let vtree_vs = T.lprp <| mkt2 lp rp Cs
  let fJs2 = T.ileaffix operator inv_op (identity 6, replicate 6 0f64) vtree_vs
  let fJs2 = map (.1) fJs2 -- Xup[i]*as'[p] + S[i]*qdd[i] + (crm vs'[i]) * vJ[i] 

  let fJs = loop fJs' = (fBs) for i < n -1 do
    let idx = n - (i+1)
    let parent = p[idx]
    in fJs' with [parent] = map2 (+) fJs'[parent] (mat_mul_vec_f64 (transpose Xup[idx]) fJs'[idx])

  in trace <| (fJs, fJs2)
  -- in trace <| map2 (\s f -> vecmul s f) S fJs  
  -- in trace test2  


def main = 
  -- let lp = sized n [0, 1, 5, 2]
  -- let rp = sized n [7, 4, 6, 3]
  let lp = [0, 1, 7, 2, 4, 8]
  let rp = [11,  6, 10, 3, 5, 9]
  let (_, p, js, _, Is, Xtrees) = autoTree 6 2 1 1
  in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3] lp rp
  -- let (_, p, js, _, Is, Xtrees) = autoTree 4 2 1 1
  -- in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, , 1] [0f64, 2, 1, 3, ] [0f64, 3, 0,  3]
  -- let (_, p, js, _, Is, Xtrees) = autoTree 100 1 1 1
  -- in rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] (replicate 100 (1f64)) (replicate 100 (1f64)) (replicate 100 (1f64))
  -- in rnea' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] [0f64, 1, 0, 0, 0, 1] [0f64, 2, 1, 3, 0, 1] [0f64, 3, 0, 0, 0, 3]


-- Ide: 1. Lav en data struktur, så man ved hvilke links er i hvilke dybder af det kinematiske træ
--      2. Kør sequentielt gennem dybderne
--      3. Ved hver dybde udregn velocity eller acceleration 
