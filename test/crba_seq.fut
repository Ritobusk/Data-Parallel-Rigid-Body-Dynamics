import "../CRBA"
import "../treeModel"
import "../lib/github.com/diku-dk/vtree/vtree"


-- Test the vtree rnea algorithm against results obtained from the matlab library made by Roy Featherstone.
-- ==
-- entry: test_crba_seq
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
entry test_crba_seq [n] (children : f64) (qs : [n]f64) (qds : [n]f64)  : ([n]f64, [n][n]f64) =
  let (_, p, js, Is, Xtrees) = autoTree n children 0 1
  let gravity = [0f64, 0, 0, 0, 0, -9.81]
  in crba_seq p js Is Xtrees gravity qs qds 
