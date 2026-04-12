import "../rnea"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"


-- Benchmark the vtree rnea algorithm.
-- ==
-- entry: test_rnea
-- input @ data/N4_bf2_rnea.in
-- output @ data/N4_bf2_rnea.out
entry test_rnea [n] (children : f64) (qs : [n]f64) (qds : [n]f64) (qdds : [n]f64)  : [n]f64 =
  let (_, p, js, Is, Xtrees) = autoTree n children 0 1
  let vtree =  T.mk_parent p (iota n)
  let tmp = T.unmk vtree
  let lp = tmp.lp 
  let rp = tmp.rp 
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in rnea'' p js Is Xtrees gravity qs qds qdds lp rp
