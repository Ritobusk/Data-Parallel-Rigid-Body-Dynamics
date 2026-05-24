
import "../rnea"
import "../rnea_optimal"
import "../treeModel"
import "../spatial_ops"
import "../lib/github.com/diku-dk/vtree/vtree"
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

entry rnea_inputC (n : i64) (bf: f64) :
    ([n]i64, [n]I_Compact, [n]X_Compact,  [n]f64, [n]f64, [n]f64) =
    let (_, p, _, Is, Xtrees) = autoTreeC n bf 0f64 1f64
    let p = sized n p
    let Is = sized n Is
    let Xtrees = sized n Xtrees

    let (q, qd, qdd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    let qdd = map f64.f32 qdd
    in (p, Is, Xtrees, q, qd, qdd )

entry rnea_input (n : i64) (bf: f64) :
    ([n]i64,  [n][6][6]f64, [n][6][6]f64, 
     [n]f64, [n]f64, [n]f64) =
    let (_, p, _, Is, Xtrees) = autoTree n bf 0f64 1f64
    let p = sized n p
    let Is = sized n Is
    let Xtrees = sized n Xtrees

    let (q, qd, qdd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    let qdd = map f64.f32 qdd
    in (p, Is, Xtrees, q, qd, qdd )

-- Benchmark the sequential with optimimal data structures rnea algorithm.
-- ==
-- entry: bench_rnea_optimized_ds_seq 
-- script input { rnea_inputC 10i64 1f64 }
-- script input { rnea_inputC 10i64 2f64 }
-- script input { rnea_inputC 10i64 10f64 }
-- script input { rnea_inputC 10i64 1000f64 }
-- script input { rnea_inputC 100i64 1f64 }
-- script input { rnea_inputC 100i64 2f64 }
-- script input { rnea_inputC 100i64 10f64 }
-- script input { rnea_inputC 100i64 1000f64 }
-- script input { rnea_inputC 1000i64 1f64 }
-- script input { rnea_inputC 1000i64 2f64 }
-- script input { rnea_inputC 1000i64 10f64 }
-- script input { rnea_inputC 1000i64 1000f64 }
-- script input { rnea_inputC 10000i64 1f64 }
-- script input { rnea_inputC 10000i64 2f64 }
-- script input { rnea_inputC 10000i64 10f64 }
-- script input { rnea_inputC 10000i64 1000f64 }
-- script input { rnea_inputC 100000i64 1f64 }
-- script input { rnea_inputC 100000i64 2f64 }
-- script input { rnea_inputC 100000i64 10f64 }
-- script input { rnea_inputC 100000i64 1000f64 }
-- script input { rnea_inputC 1000000i64 1f64 }
-- script input { rnea_inputC 1000000i64 2f64 }
-- script input { rnea_inputC 1000000i64 10f64 }
-- script input { rnea_inputC 1000000i64 1000f64 }
-- script input { rnea_inputC 2000000i64 1f64 }
-- script input { rnea_inputC 4000000i64 1f64 }
entry bench_rnea_optimized_ds_seq [n] (p : [n]i64)  (Is : [n]I_Compact) (Xtrees: [n]X_Compact) (q : [n]f64)  (qd : [n]f64) (qdd : [n]f64) : [n]f64 =
  let gravity = {w = [0,0,0f64], v_O = [0,0, -9.81f64]}
  in rnea_optimized_ds_seq p (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd qdd

-- Benchmark the sequential rnea algorithm.
-- ==
-- entry: bench_rnea_no_vtree
-- script input { rnea_input 10i64 1f64 }
-- script input { rnea_input 10i64 2f64 }
-- script input { rnea_input 10i64 10f64 }
-- script input { rnea_input 10i64 1000f64 }
-- script input { rnea_input 100i64 1f64 }
-- script input { rnea_input 100i64 2f64 }
-- script input { rnea_input 100i64 10f64 }
-- script input { rnea_input 100i64 1000f64 }
-- script input { rnea_input 1000i64 1f64 }
-- script input { rnea_input 1000i64 2f64 }
-- script input { rnea_input 1000i64 10f64 }
-- script input { rnea_input 1000i64 1000f64 }
-- script input { rnea_input 10000i64 1f64 }
-- script input { rnea_input 10000i64 2f64 }
-- script input { rnea_input 10000i64 10f64 }
-- script input { rnea_input 10000i64 1000f64 }
-- script input { rnea_input 100000i64 1f64 }
-- script input { rnea_input 100000i64 2f64 }
-- script input { rnea_input 100000i64 10f64 }
-- script input { rnea_input 100000i64 1000f64 }
-- script input { rnea_input 1000000i64 1f64 }
-- script input { rnea_input 1000000i64 2f64 }
-- script input { rnea_input 1000000i64 10f64 }
-- script input { rnea_input 1000000i64 1000f64 }
-- script input { rnea_input 2000000i64 1f64 }
-- script input { rnea_input 4000000i64 1f64 }
entry bench_rnea_no_vtree [n] (p : [n]i64)  (Is : [n][6][6]f64) (Xtrees: [n][6][6]f64) (q : [n]f64)  (qd : [n]f64) (qdd : [n]f64) : [n]f64 =
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in rnea p (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd qdd
