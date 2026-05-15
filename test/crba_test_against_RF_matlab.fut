import "../CRBA"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"


-- Test the vtree rnea algorithm against results obtained from the matlab library made by Roy Featherstone.
-- ==
-- entry: test_crba
-- input @ data/N5_bf1_crba.in
-- output @ data/N5_bf1_crba.out
-- input @ data/N5_bf2_crba.in
-- output @ data/N5_bf2_crba.out
-- input @ data/N500_bf1_crba.in
-- output @ data/N500_bf1_crba.out
-- input @ data/N500_bf100_crba.in
-- output @ data/N500_bf100_crba.out
-- input @ data/N1000_bf1000_crba.in
-- output @ data/N1000_bf1000_crba.out
entry test_crba [n] (children : f64) (qs : [n]f64) (qds : [n]f64)  : ([n]f64, [n][n]f64) =
  let (_, _, js, Is, Xtrees, lp, rp, paths, p_ii1) = autoVTree n children 0 1
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in crba_vtree_optimized js Is Xtrees gravity qs qds lp rp paths p_ii1
