import "../rnea"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"


-- Test the vtree rnea algorithm against results obtained from the matlab library made by Roy Featherstone.
-- ==
-- entry: test_rnea
-- input @ data/N4_bf2_rnea.in
-- output @ data/N4_bf2_rnea.out
-- input @ data/N6_bf2_rnea.in
-- output @ data/N6_bf2_rnea.out
-- input @ data/N10_bf1_rnea.in
-- output @ data/N10_bf1_rnea.out
-- input @ data/N100_bf2_rnea.in
-- output @ data/N100_bf2_rnea.out
entry test_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64)  : [n]f64 =
  let (_, p, js, _, Is, Xtrees) = autoTree n children 1 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  in rnea' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] qs qds qdds 
  -- in and <| map2 (\x y -> f64.abs (x - y) < 1e-8f64) vtree_res matlab_res

