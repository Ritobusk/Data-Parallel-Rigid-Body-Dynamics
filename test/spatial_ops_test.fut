import "../treeModel"
import "../spatial_ops"
import "../matrix_ops"
import "../lib/github.com/diku-dk/cpprandom/random"

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

-- Tests the optimized/compact data structures against the normal 6d vectors and 6x6 matrices of spatial vector algebra


def transform_multiplication (n : i64) : bool =
  let (q, _, _) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
  let q = map f64.f32 q

  let (_, p, js, _, Xtrees) = autoTree n 1 1 1
  let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js q 
  let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtrees

  let (_, _, jsC, _, XtreesC) = autoTreeC n 1 1 1
  let (XJC, SC) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) jsC q 
  let XupC      = map2 (\xj xtree -> transform_XX xj xtree) XJC XtreesC
  let XupC' = map (from_XC_to_aXb) XupC

  in map2 (\xup xup' -> map2 (\r r' -> map2 (\x y -> f64.abs (x - y) < error_tolerance) r r') xup xup') Xup XupC'
      |> flatten |> flatten
      |> reduce (&&) true


-- ==
-- entry: plucker_transform 
-- input {4i64} output {true}
-- input {8i64} output {true}
entry plucker_transform (n : i64) =
  transform_multiplication n
