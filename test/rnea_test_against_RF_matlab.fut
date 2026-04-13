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
-- input @ data/N1000_bf2_rnea.in
-- output @ data/N1000_bf2_rnea.out
-- input @ data/N1000_bf100_rnea.in
-- output @ data/N1000_bf100_rnea.out
entry test_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64)  : [n]f64 =
  let (_, p, js, Is, Xtrees) = autoTree n children 0 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in rnea'' js Is Xtrees gravity qs qds qdds lp rp



-- Test the vtree rnea algorithm with external forces.
-- ==
-- entry: test_fext_rnea
-- input  @ data/10N_2bf_fext_rnea.in
-- output @ data/10N_2bf_fext_rnea.out
-- input  @ data/10N_1bf_fext_rnea.in
-- output @ data/10N_1bf_fext_rnea.out
-- input  @ data/50N_7bf_fext_rnea.in
-- output @ data/50N_7bf_fext_rnea.out
entry test_fext_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64) (fext : [n][6]f64)  : [n]f64 =
  let (_, p, js, Is, Xtrees) = autoTree n children 0 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in rnea_vtree_with_f_ext js Is Xtrees gravity qs qds qdds fext lp rp
