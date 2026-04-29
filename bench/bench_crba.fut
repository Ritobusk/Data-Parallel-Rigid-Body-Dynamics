import "../CRBA"
import "../treeModel"
import "../spatial_ops"
import "../lib/github.com/diku-dk/vtree/vtree"

entry rnea_input [n] (a : (i64, f64, f64, f64, [n]f64, [n]f64) ) :   
    ([n][6][6]f64, [n][6][6]f64, [n]i64, [n]i64, []i64, []i64,
     [n]f64, [n]f64) =
    let (_, _, _, Is, Xtrees, lp, rp, paths, p_ii1) = autoVTree a.0 a.1 a.2 a.3
    let Is = sized n Is
    let Xtrees = sized n Xtrees
    let lp = sized n lp
    let rp = sized n rp
    in (Is, Xtrees, lp, rp, paths, p_ii1, a.4, a.5)


-- Benchmark the vtree rnea algorithm.
-- ==
-- entry: bench_crba
-- script input { rnea_input ($loaddata "data/N10_bf1_crba.in") }
-- script input { rnea_input ($loaddata "data/N10_bf2_crba.in") }
-- script input { rnea_input ($loaddata "data/N10_bf100_crba.in") }
-- script input { rnea_input ($loaddata "data/N50_bf1_crba.in") }
-- script input { rnea_input ($loaddata "data/N50_bf2_crba.in") }
-- script input { rnea_input ($loaddata "data/N50_bf100_crba.in") }
-- script input { rnea_input ($loaddata "data/N100_bf1_crba.in") }
-- script input { rnea_input ($loaddata "data/N100_bf2_crba.in") }
-- script input { rnea_input ($loaddata "data/N100_bf100_crba.in") }
-- script input { rnea_input ($loaddata "data/N500_bf1_crba.in") }
-- script input { rnea_input ($loaddata "data/N500_bf2_crba.in") }
-- script input { rnea_input ($loaddata "data/N500_bf100_crba.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf1_crba.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf2_crba.in") }
-- script input { rnea_input ($loaddata "data/N5000_bf100_crba.in") }
entry bench_crba [n] [nd]  (Is : [n][6][6]f64) (Xtrees: [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) (paths : [nd]i64) (p_ii1 : [nd]i64) (q : [n]f64)  (qd : [n]f64)  : ([n]f64, [n][n]f64) =
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in crba_vtree' (replicate n #Rz : [n]jointT) Is Xtrees gravity q qd lp rp paths p_ii1

