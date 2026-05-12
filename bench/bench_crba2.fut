import "../CRBA"
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
    let (_, qds) = unzip (map (dist.rand d) rngs')
    in (qs, qds) 
}

module test_f32_rand_m =
  mktest (uniform_real_distribution f32 minstd_rand)

entry crba_input (n : i64) (bf: f64) :
    ([n][6][6]f64, [n][6][6]f64, [n]i64, [n]i64, []i64, []i64,
     [n]f64, [n]f64) =
    let (_, _, _, Is, Xtrees, lp, rp, paths, p_ii1) = autoVTree n bf 0f64 1f64
    let Is = sized n Is
    let Xtrees = sized n Xtrees
    let lp = sized n lp
    let rp = sized n rp
    let (q, qd) =  test_f32_rand_m.test (123i32) n (-1f32, 1f32)
    let q = map f64.f32 q
    let qd = map f64.f32 qd
    in (Is, Xtrees, lp, rp, paths, p_ii1, q, qd)


-- Benchmark the vtree rnea algorithm.
-- ==
-- entry: bench_crba
-- script input { crba_input 10i64 1f64 }
-- script input { crba_input 10i64 2f64 }
-- script input { crba_input 10i64 10f64 }
-- script input { crba_input 10i64 1000f64 }
-- script input { crba_input 100i64 1f64 }
-- script input { crba_input 100i64 2f64 }
-- script input { crba_input 100i64 10f64 }
-- script input { crba_input 100i64 1000f64 }
-- script input { crba_input 1000i64 1f64 }
-- script input { crba_input 1000i64 2f64 }
-- script input { crba_input 1000i64 10f64 }
-- script input { crba_input 1000i64 1000f64 }
-- script input { crba_input 5000i64 1f64 }
-- script input { crba_input 5000i64 2f64 }
-- script input { crba_input 5000i64 10f64 }
-- script input { crba_input 5000i64 1000f64 }
entry bench_crba [n] [nd]  (Is : [n][6][6]f64) (Xtrees: [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) (paths : [nd]i64) (p_ii1 : [nd]i64) (q : [n]f64)  (qd : [n]f64)  : ([n]f64, [n][n]f64) =
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in crba_vtree' (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd lp rp paths p_ii1

