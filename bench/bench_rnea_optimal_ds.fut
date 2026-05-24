
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

entry rnea_input (n : i64) (bf: f64) :
    ([n]I_Compact, [n]X_Compact, [n]i64, [n]i64, [n]f64, [n]f64, [n]f64) =
    let (_, p, _, Is, Xtrees) = autoTreeC n bf 0f64 1f64
    let p = sized n p
    let Is = sized n Is
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 

    let (q, qd, qdd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    let qdd = map f64.f32 qdd
    in (Is, Xtrees, lp, rp, q, qd, qdd )


-- Benchmark the vtree rnea algorithm.
-- ==
-- entry: bench_rnea
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
entry bench_rnea [n] (Is : [n]I_Compact) (Xtrees: [n]X_Compact) (lp : [n]i64) (rp : [n]i64) (q : [n]f64)  (qd : [n]f64) (qdd : [n]f64) : [n]f64 =
  let gravity = {w = [0,0,0f64], v_O = [0,0, -9.81f64]}
  in rnea_vtree_optimized_ds (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd qdd lp rp

