import "../rnea_optimal"
import "../spatial_ops"
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
-- input @ data/N1000_bf2_rnea.in
-- output @ data/N1000_bf2_rnea.out
-- input @ data/N1000_bf100_rnea.in
-- output @ data/N1000_bf100_rnea.out
-- input @ data/N100000_bf2_rnea.in
-- output @ data/N100000_bf2_rnea.out
entry test_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64)  : [n]f64 =
  let (_, p, js, Is, Xtrees) = autoTreeC n children 0 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  let gravity = {w = [0,0,0f64], v_O = [0,0, -9.81f64]}
  in rnea_vtree_optimized_ds js Is Xtrees gravity qs qds qdds lp rp
