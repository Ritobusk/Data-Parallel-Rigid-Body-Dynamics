import "../rnea"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"


-- Test the vtree rnea algorithm against results obtained from the matlab library made by Roy Featherstone.
-- ==
-- entry: test_rnea
-- input @ data/N4_bf2_rnea.in
-- output { true }
-- input @ data/N6_bf2_rnea.in
-- output { true }
entry test_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64) (matlab_res : [n]f64) : bool =
  let (_, p, js, _, Is, Xtrees) = autoTree n children 1 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  let vtree_res = rnea'' p js Is Xtrees [0f64, 0, 0, 0, 0, -9.81] qs qds qdds lp rp
  in and <| map2 (\x y -> f64.abs (x - y) < 1e-8f64) vtree_res matlab_res

