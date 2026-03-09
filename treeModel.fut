import "matrix_ops"
import "spatial_ops"


-- Inspiration is taken from: https://royfeatherstone.org/spatial/v2/sourceText/autoTree.txt
-- Creates a kinematic tree 
def autoTree (nb : i64) (bf : f64) (skew : f64) (taper : f64) =
    let ids = iota nb
    let joint_types = replicate nb (#Rx) : [nb]jointT
    let parents = map (\i -> i64.f64 <| (f64.floor ( (((f64.i64 i) + 1.0) - 2.0 + (f64.ceil bf) ) / bf )) - 1.0 ) ids

    let lengths = map (\i -> taper ** (f64.i64 i)) ids 
    let masses = map (\i ->  taper ** (3 * (f64.i64 i))) ids 
    let CoMs   = map (\l -> [0.5 * l, 0, 0]) lengths
    let Icms   = map (\i ->
                  let m = masses[i]
                  let l = lengths[i]
                  let d = diagonal [0.0025,1.015/12,1.015/12]  -- each body is a cylinder
                  in scal_mul_mat (m * l) d) ids

    let Is = map (\i -> mcI masses[i] CoMs[i] Icms[i]) ids 

    let Xtree       = map (\i -> if i == 0 then xlt [0,0,0] -- identity
                                 else matmul_f64 (rotx skew)  (xlt [lengths[parents[i]], 0, 0]) -- Not sure if this just puts all the children in the excact same place or not...
                          ) ids

    in (ids, parents, joint_types,  lengths, Is, Xtree)

def main =
  autoTree 2 1 1 1
