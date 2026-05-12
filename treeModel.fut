import "matrix_ops"
import "spatial_ops"
import "lib/github.com/diku-dk/vtree/vtree"
import "lib/github.com/diku-dk/segmented/segmented"

module T = vtree

-- Inspiration is taken from: https://royfeatherstone.org/spatial/v2/sourceText/autoTree.txt
-- Creates a kinematic tree 
def autoTree (nb : i64) (bf : f64) (skew : f64) (taper : f64) :
    ([nb]i64, [nb]i64, [nb]jointT, [nb][6][6]f64, [nb][6][6]f64) =
    let ids = iota nb
    let joint_types = replicate nb (#Rz) : [nb]jointT
    let parents = map (\i -> i64.f64 <| (f64.floor ( (((f64.i64 i) + 1.0) - 2.0 + (f64.ceil bf) ) / bf )) - 1.0 ) ids
    let parents' = parents with [0] = 0  -- The root's parent is itself.

    let lengths = map (\i -> taper ** (f64.i64 i)) ids 
    let masses = map (\i ->  taper ** (3 * (f64.i64 i))) ids 
    let CoMs   = map (\l -> [0.5 * l, 0, 0]) lengths
    let Icms   = map (\i ->
                  let m = masses[i]
                  let l = lengths[i]
                  let d = diagonal [0.0025,1.015/12,1.015/12]  -- each body is a cylinder
                  in scal_mul_mat_f64 (m * (l**2)) d) ids

    let Is = map (\i -> mcI masses[i] CoMs[i] Icms[i]) ids 

    let Xtree       = map (\i -> if i == 0 then xlt [0,0,0] -- identity
                                 else matmul_f64 (rotx skew)  (xlt [lengths[parents'[i]], 0, 0]) -- Not sure if this just puts all the children in the excact same place or not...
                          ) ids

    in (ids, parents', joint_types, Is, Xtree)

-- Creates model data with the lp and rp for a vtree
def autoVTree (nb : i64) (bf : f64) (skew : f64) (taper : f64) :
    ([nb]i64, [nb]i64, [nb]jointT, [nb][6][6]f64, [nb][6][6]f64, [nb]i64, [nb]i64, []i64, []i64) =
    let ids = iota nb
    let joint_types = replicate nb (#Rz) : [nb]jointT
    let parents = map (\i -> i64.f64 <| (f64.floor ( (((f64.i64 i) + 1.0) - 2.0 + (f64.ceil bf) ) / bf )) - 1.0 ) ids
    let parents' = parents with [0] = 0  -- The root's parent is itself.

    let lengths = map (\i -> taper ** (f64.i64 i)) ids 
    let masses = map (\i ->  taper ** (3 * (f64.i64 i))) ids 
    let CoMs   = map (\l -> [0.5 * l, 0, 0]) lengths
    let Icms   = map (\i ->
                  let m = masses[i]
                  let l = lengths[i]
                  let d = diagonal [0.0025,1.015/12,1.015/12]  -- each body is a cylinder
                  in scal_mul_mat_f64 (m * (l**2)) d) ids

    let Is = map (\i -> mcI masses[i] CoMs[i] Icms[i]) ids 

    let Xtree       = map (\i -> if i == 0 then xlt [0,0,0] -- identity
                                 else matmul_f64 (rotx skew)  (xlt [lengths[parents'[i]], 0, 0]) -- Not sure if this just puts all the children in the excact same place or not...
                          ) ids

    let vtree = T.mk_parent parents' (iota nb)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 

    let depths = T.depth vtree
    let tree_depth = reduce i64.max 0i64 depths
    let p_ii1 = replicated_iota depths 
    let (_, paths) = map2 (\p d -> 
        loop (j, path) = (p, replicate tree_depth (-1i64) ) for k < (d) do
                let path = path with [k] = j 
                let j = parents'[j]
                in (j, path)
        ) parents' depths 
        |> unzip
    let paths = flatten paths |> filter (> -1i64) 

    in (ids, parents', joint_types, Is, Xtree, lp, rp, paths, p_ii1)
