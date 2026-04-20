import "../rnea"
import "../treeModel"
import "../spatial_ops"
import "../lib/github.com/diku-dk/vtree/vtree"

entry rnea_input [n] (a : (i64, f64, f64, f64, [n]f64, [n]f64, [n]f64) ) :   
    ([n][6][6]f64, [n][6][6]f64, [n]i64, [n]i64,
     [n]f64, [n]f64, [n]f64) =
    let (_, p, _, Is, Xtrees) = autoTree a.0 a.1 a.2 a.3
    let p = sized n p
    let Is = sized n Is
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    in (Is, Xtrees, lp, rp, a.4, a.5, a.6)


-- Benchmark the vtree rnea algorithm.
-- ==
-- entry: bench_rnea
-- script input { rnea_input ($loaddata "data/test.in") }
-- script input { rnea_input ($loaddata "data/N100_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N10000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N50000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N100000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000000_bf1_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N100_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N10000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N50000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N100000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000000_bf2_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N100_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N10000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N50000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N100000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N500000_bf100_rnea_bench.in") }
-- script input { rnea_input ($loaddata "data/N1000000_bf100_rnea_bench.in") }
entry bench_rnea [n]  (Is : [n][6][6]f64) (Xtrees: [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) (q : [n]f64)  (qd : [n]f64) (qdd : [n]f64) : [n]f64 =
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in rnea_vtree_optimized (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd qdd lp rp

