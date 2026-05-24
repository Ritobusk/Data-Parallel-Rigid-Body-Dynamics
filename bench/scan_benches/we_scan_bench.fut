import "../../matrix_ops"
import "../../spatial_ops"
import "../../treeModel"
import "scan_variations"
import "vtree_with_scans"

module T = vtree

def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

entry vtree_vectoradd  (n : i64) :   
    ([n][6]f64, [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTree n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js (replicate n 1f64) 
    let vJ      = map2 (\s v -> map (\x -> x * v) s) S (replicate n 1f64) 
    in (vJ, lp, rp)

entry vtree_matrixmul  (n : i64) :   
    ([n][6][6]f64, [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTree n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js (replicate n 1f64) 
    let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtrees
    in (Xup, lp, rp)

entry vtree_matrixmul_C  (n : i64) :   
    ([n]X_Compact, [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTreeC n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) js (replicate n 1f64) 
    let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtrees
    in (Xup, lp, rp)

entry complex_scan_input  (n : i64) :   
    ([n]([6][6]f64, [6]f64), [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTree n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalc joint j_pos) js (replicate n 1f64) 
    let vJ      = map2 (\s v -> map (\x -> x * v) s) S (replicate n 1f64) 
    let Xup     = map2 (\xj xtree -> matmul_f64 xj xtree) XJ Xtrees
    let Cs = zip Xup vJ
    in (Cs, lp, rp)

entry complex_scan_input_C  (n : i64) :   
    ([n](X_Compact, mv), [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTreeC n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) js (replicate n 1f64) 
    let vJ      = map2 (\s v ->  s `scal_mul_mv` v) (replicate n 1f64) S 
    let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtrees
    let Cs = zip Xup vJ
    in (Cs, lp, rp)

-- ==
-- entry: bench_we_rootfix_va
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
entry bench_we_rootfix_va [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_we (vecadd_f64) (scal_mul_vec_f64 (-1)) (replicate 6 0f64) t

-- ==
-- entry: bench_we_rootfix_mm
-- script input { vtree_matrixmul  100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
entry bench_we_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_we (matmul_f64) (XBtoA_from_XAtoB_M)  (identity 6) t

-- ==
-- entry: bench_we_rootfix_mm_C
-- script input { vtree_matrixmul_C  100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
-- script input { vtree_matrixmul_C 4000000i64 }  
-- script input { vtree_matrixmul_C 8000000i64 }  
entry bench_we_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_we (transform_XX) (transform_inv) (copy transform_identity) t 

-- ==
-- entry: test_work_efficient_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
entry test_work_efficient_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_we operator inv_op (identity 6, replicate 6 0f64) t

-- ==
-- entry: test_work_efficient_rootfix_C
-- script input { complex_scan_input_C 100i64 }  
-- script input { complex_scan_input_C 1000i64 }  
-- script input { complex_scan_input_C 10000i64 }  
-- script input { complex_scan_input_C 100000i64 }  
-- script input { complex_scan_input_C 1000000i64 }  
-- script input { complex_scan_input_C 2000000i64 }  
-- script input { complex_scan_input_C 4000000i64 }  
-- script input { complex_scan_input_C 8000000i64 }  
entry test_work_efficient_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_we operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 
