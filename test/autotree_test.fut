import "../treeModel"

-- Test the auto tree against Roy Featherstone's auto tree.
-- ==
-- entry: test_autoTree
-- input  @ data/10_1_at.in
-- output @ data/10_1_at.out
-- input  @ data/10_4_at.in
-- output @ data/10_4_at.out
-- input  @ data/100_2_at.in
-- output @ data/100_2_at.out
-- input  @ data/100_1point5_at.in
-- output @ data/100_1point5_at.out
entry test_autoTree  (n : i64) (children : f64) (skew : f64) (taper : f64)  
  : ([n]i64, [n][6][6]f64, [n][6][6]f64)  =
  let (_, p, _, Is, Xtrees) = autoTree n children skew taper
  in (map (\i -> if i == 0 then p[i] else p[i]+1) (indices p),
      Xtrees,
      Is)
