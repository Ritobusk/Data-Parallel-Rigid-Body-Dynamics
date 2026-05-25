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
-- entry: test_blocked_512_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
entry test_blocked_512_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (transform_XX) (transform_inv) (copy transform_identity) t 512

-- ==
-- entry: test_blocked_256_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
entry test_blocked_256_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (transform_XX) (transform_inv) (copy transform_identity) t 256

-- ==
-- entry: test_blocked_128_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
entry test_blocked_128_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (transform_XX) (transform_inv) (copy transform_identity) t 128

-- ==
-- entry: test_blocked_64_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
entry test_blocked_64_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (transform_XX) (transform_inv) (copy transform_identity) t 64

-- ==
-- entry: test_blocked_32_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
entry test_blocked_32_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (transform_XX) (transform_inv) (copy transform_identity) t 32


-- ==
-- entry: test_blocked_512_rootfix_C
-- script input { complex_scan_input_C  100i64 }  
-- script input { complex_scan_input_C  1000i64 }  
-- script input { complex_scan_input_C  10000i64 }  
-- script input { complex_scan_input_C  100000i64 }  
-- script input { complex_scan_input_C  1000000i64 }  
-- script input { complex_scan_input_C  2000000i64 }  
entry test_blocked_512_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 512

-- ==
-- entry: test_blocked_256_rootfix_C
-- script input { complex_scan_input_C  100i64 }  
-- script input { complex_scan_input_C  1000i64 }  
-- script input { complex_scan_input_C  10000i64 }  
-- script input { complex_scan_input_C  100000i64 }  
-- script input { complex_scan_input_C  1000000i64 }  
-- script input { complex_scan_input_C  2000000i64 }  
entry test_blocked_256_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 256

-- ==
-- entry: test_blocked_128_rootfix_C
-- script input { complex_scan_input_C  100i64 }  
-- script input { complex_scan_input_C  1000i64 }  
-- script input { complex_scan_input_C  10000i64 }  
-- script input { complex_scan_input_C  100000i64 }  
-- script input { complex_scan_input_C  1000000i64 }  
-- script input { complex_scan_input_C  2000000i64 }  
entry test_blocked_128_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 128

-- ==
-- entry: test_blocked_64_rootfix_C
-- script input { complex_scan_input_C  100i64 }  
-- script input { complex_scan_input_C  1000i64 }  
-- script input { complex_scan_input_C  10000i64 }  
-- script input { complex_scan_input_C  100000i64 }  
-- script input { complex_scan_input_C  1000000i64 }  
-- script input { complex_scan_input_C  2000000i64 }  
entry test_blocked_64_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 64

-- ==
-- entry: test_blocked_32_rootfix_C
-- script input { complex_scan_input_C  100i64 }  
-- script input { complex_scan_input_C  1000i64 }  
-- script input { complex_scan_input_C  10000i64 }  
-- script input { complex_scan_input_C  100000i64 }  
-- script input { complex_scan_input_C  1000000i64 }  
-- script input { complex_scan_input_C  2000000i64 }  
entry test_blocked_32_rootfix_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact,mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t 32

-- ==
-- entry: test_blocked_128_rootfix_va
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
-- script input { vtree_vectoradd 2000000i64 }  
entry test_blocked_128_rootfix_va [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (vecadd_f64) (scal_mul_vec_f64 (-1)) (replicate 6 0f64) t 128

-- ==
-- entry: test_blocked_64_rootfix_va
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
-- script input { vtree_vectoradd 2000000i64 }  
entry test_blocked_64_rootfix_va [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (vecadd_f64) (scal_mul_vec_f64 (-1)) (replicate 6 0f64) t 64

-- ==
-- entry: test_blocked_512_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
entry test_blocked_512_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t 512
-- ==
-- entry: test_blocked_256_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
entry test_blocked_256_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t 256

-- ==
-- entry: test_blocked_128_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
entry test_blocked_128_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t 128


-- ==
-- entry: test_blocked_64_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
entry test_blocked_64_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t 64

-- ==
-- entry: test_blocked_32_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
entry test_blocked_32_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t 32


-- ==
-- entry: test_blocked_512_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_512_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 512i64

-- ==
-- entry: test_blocked_256_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_256_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 256i64

-- ==
-- entry: test_blocked_128_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_128_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 128i64

-- ==
-- entry: test_blocked_64_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_64_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 64i64

-- ==
-- entry: test_blocked_32_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_32_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 32i64

-- ==
-- entry: test_blocked_16_rootfix
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
entry test_blocked_16_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 16i64


