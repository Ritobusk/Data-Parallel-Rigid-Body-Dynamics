import "../../matrix_ops"
import "../../spatial_ops"
import "../../treeModel"
import "vtree_with_scans"

module T = vtree

def mkt2 'a [n] (lp: [n]i64) (rp: [n]i64) (ds:[n]a) : {lp: [n]i64, rp: [n]i64, data: [n]a} =
    {lp=lp,rp=rp,data=ds}

-- Blocked scan from Troels
def blocked 'a [n] (op: a -> a -> a) (ne: a) (xs: [n]a) (bs : i64) : [n]a =
    let block_size = #[param(scan_block_size)] bs
    let num_blocks = (n + block_size - 1) / block_size
    let block_scans =
      tabulate num_blocks (\i ->
        let block =
          tabulate block_size (\j ->
            let l = i * block_size + j
            in if l < n then xs[l] else ne)
        in #[sequential] scan op ne block)
    let carry_outs =
      -- Can use any scan for this.
      scan op ne (map last block_scans)
    in map2 (\x l ->
               let i = l / block_size
               in if i > 0 then carry_outs[i - 1] `op` x else x)
          (take n (flatten block_scans))
          (iota n)

def exscan_blocked f ne xs bs =
  map2 (\i x -> if i == 0 then ne else x)
       (indices xs)
       (rotate (-1) (blocked f ne xs bs))

-- The scan implementations are my DPP implementation where I modified them to take
-- an operator as input
def ilog2 (x: i64 ) = 63 - i64.i32 (i64.clz x)

def next_power_of_two (n : i64) : i64 =
  loop acc = 2 while acc < n do 
    acc * 2

def exscan f ne xs =
    map2 (\i x -> if i == 0 then ne else x)
         (indices xs)
         (rotate (-1) (scan  f ne xs))

def work_efficient [n] 'a (op : a -> a -> a) (ne : a)  (xs : [n]a)  : [n]a =
    let k = next_power_of_two n 
    let m = ilog2 k
    let xs' = tabulate k (\i -> if i < n then xs[i] else ne)
    let upswept =
        loop xs' = copy xs' for d in 0...m-1 do
            let offset = 2 ** (d + 1)
            let indx_to_update = iota (k/(offset)) 
                                 |> map (\x -> (x + 1) * offset - 1)
            let update_vals = map (\i -> op xs'[i - (offset / 2)] xs'[i] ) indx_to_update
            in (scatter xs' indx_to_update update_vals)
    let upswept[k-1] = ne
    let downswept =
        loop xs' = copy (upswept) for d in m..m-1...1 do
            let offset = 2 ** (d - 1)
            let indx_to_update = iota (k/(offset)) 
                                    |> map (\x -> (x + 1) * offset - 1)
            let update_vals = 
                map2 
                    (\itu i -> if i % 2 == 1 then op xs'[itu] xs'[itu - (offset)]
                                             else xs'[itu + offset])
                    indx_to_update (indices indx_to_update)
            in (scatter xs' indx_to_update update_vals)
    in take n downswept

def inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
  let inv_cia = XBtoA_from_XAtoB_M ci.0
  let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
  in (inv_cia, inv_cib)

def operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
  (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)

def inv_op_C (ci : (X_Compact, mv)) : (X_Compact, mv) =
  let inv_cia = transform_inv ci.0
  let inv_cib = scal_mul_mv (-1) (inv_cia `Xm` ci.1)
  in (inv_cia, inv_cib)

def operator_C (si : (X_Compact, mv)) (ci : (X_Compact, mv)) : (X_Compact, mv) =
  (ci.0 `transform_XX` si.0,    (ci.0 `Xm` si.1) `mv_add` ci.1)

def rootfix_work_efficient 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = work_efficient op ne R 
    in map (\i -> S[i]) lp

def rootfix_vector_add [n]   (lp : [n]i64) (rp : [n]i64) (data : [n][6]f64) : [n][6]f64 =
    let I = replicate (2 * n) (replicate 6 0f64)
    let L = scatter I lp data
    let R = scatter L rp (map ((scal_mul_vec_f64 (-1))) data)
    let S0 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][0])) 
    let S1 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][1])) 
    let S2 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][2])) 
    let S3 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][3])) 
    let S4 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][4])) 
    let S5 = scan (+) 0 (iota (2*n) |> map (\i -> R[i][5])) 
    in map (\i -> [S0[i], S1[i], S2[i], S3[i], S4[i], S5[i] ]) lp

def rootfix_vector_add2 [n]   (lp : [n]i64) (rp : [n]i64) (data : [n][6]f64) : [n][6]f64 =
    let I = replicate (2 * n) (replicate 6 0f64)
    let L = scatter I lp data
    let R = scatter L rp (map ((scal_mul_vec_f64 (-1))) data)
    let S = transpose R
            |> map (scan (+) 0f64)
            |> transpose
    in map (\i -> S[i]) lp

def irootfix_vector_add2 'a [n] (lp: [n]i64) (rp: [n]i64) (data : [n][6]f64) : [n][6]f64 =
    map2 vecadd_f64 (rootfix_vector_add2 lp rp data) data

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

entry vtree_vectoraddC  (n : i64) :   
    ([n]mv, [n]X_Compact, [n]i64, [n]i64) =
    let (_, p, js, _, Xtrees) = autoTreeC n 2 0 1
    let p = sized n p
    let Xtrees = sized n Xtrees
    let vtree =  T.mk_parent p  (iota n)
    let tmp = T.unmk vtree
    let lp = tmp.lp 
    let rp = tmp.rp 
    let (XJ, S) = unzip <| map2 (\joint j_pos -> jcalcC joint j_pos) js (replicate n 1f64) 
    let vJ      = map2 (scal_mul_mv) (replicate n 1f64) S
    let Xup     = map2 (\xj xtree -> transform_XX xj xtree) XJ Xtrees
    in (vJ, Xup, lp, rp)


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
-- entry: test_unfolded_rootfix_va_with_conversion
-- script input { vtree_vectoraddC 100i64 }  
-- script input { vtree_vectoraddC 1000i64 }  
-- script input { vtree_vectoraddC 10000i64 }  
-- script input { vtree_vectoraddC 100000i64 }  
-- script input { vtree_vectoraddC 1000000i64 }  
-- script input { vtree_vectoraddC 2000000i64 }  
-- script input { vtree_vectoraddC 4000000i64 }  
-- script input { vtree_vectoraddC 8000000i64 }  
entry test_unfolded_rootfix_va_with_conversion [n] (vJ : [n]mv) (Xup : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]mv=
    let vtree_transform = T.lprp <| mkt2 lp rp Xup
    let transformation_tree = T.irootfix_blocked transform_XX_rev X_inv (copy transform_identity) vtree_transform 64i64
    let vJ'  = map2 Xm transformation_tree vJ
    let vJ'' = map mv_to_6d vJ'
    in irootfix_vector_add2  lp rp vJ''
        |> map d6_to_mv
        

-- ==
-- entry: test_normal_leaffix_va
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
-- script input { vtree_vectoradd 2000000i64 }  
-- script input { vtree_vectoradd 4000000i64 }  
-- script input { vtree_vectoradd 8000000i64 }  
entry test_normal_leaffix_va [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.ileaffix (vecadd_f64) (scal_mul_vec_f64 (-1)) (replicate 6 0f64) t

-- ==
-- entry: test_unfolded_rootfix_va2
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
-- script input { vtree_vectoradd 2000000i64 }  
-- script input { vtree_vectoradd 4000000i64 }  
-- script input { vtree_vectoradd 8000000i64 }  
entry test_unfolded_rootfix_va2 [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
     rootfix_vector_add2  lp rp data

-- ==
-- entry: test_unfolded_rootfix_va
-- script input { vtree_vectoradd 100i64 }  
-- script input { vtree_vectoradd 1000i64 }  
-- script input { vtree_vectoradd 10000i64 }  
-- script input { vtree_vectoradd 100000i64 }  
-- script input { vtree_vectoradd 1000000i64 }  
-- script input { vtree_vectoradd 2000000i64 }  
-- script input { vtree_vectoradd 4000000i64 }  
-- script input { vtree_vectoradd 8000000i64 }  
entry test_unfolded_rootfix_va [n] (data : [n][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6]f64 =
     rootfix_vector_add  lp rp data

-- ==
-- entry: test_normal_rootfix_mm
-- script input { vtree_matrixmul 100i64 }  
-- script input { vtree_matrixmul 1000i64 }  
-- script input { vtree_matrixmul 10000i64 }  
-- script input { vtree_matrixmul 100000i64 }  
-- script input { vtree_matrixmul 1000000i64 }  
-- script input { vtree_matrixmul 2000000i64 }  
-- script input { vtree_matrixmul 4000000i64 }  
-- script input { vtree_matrixmul 8000000i64 }  
entry test_normal_rootfix_mm [n] (data : [n][6][6]f64) (lp : [n]i64) (rp : [n]i64) : [n][6][6]f64 =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix (matmul_f64) (XBtoA_from_XAtoB_M ) (identity 6) t

-- ==
-- entry: test_normal_rootfix_mm_C
-- script input { vtree_matrixmul_C 100i64 }  
-- script input { vtree_matrixmul_C 1000i64 }  
-- script input { vtree_matrixmul_C 10000i64 }  
-- script input { vtree_matrixmul_C 100000i64 }  
-- script input { vtree_matrixmul_C 1000000i64 }  
-- script input { vtree_matrixmul_C 2000000i64 }  
-- script input { vtree_matrixmul_C 4000000i64 }  
-- script input { vtree_matrixmul_C 8000000i64 }  
entry test_normal_rootfix_mm_C [n] (data : [n]X_Compact) (lp : [n]i64) (rp : [n]i64) : [n]X_Compact =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix (transform_XX) (transform_inv) (copy transform_identity) t

-- ==
-- entry: test_rootfix_normal_scan
-- script input { complex_scan_input 100i64 }  
-- script input { complex_scan_input 1000i64 }  
-- script input { complex_scan_input 10000i64 }  
-- script input { complex_scan_input 100000i64 }  
-- script input { complex_scan_input 1000000i64 }  
-- script input { complex_scan_input 2000000i64 }  
-- script input { complex_scan_input 4000000i64 }  
-- script input { complex_scan_input 8000000i64 }  
entry test_rootfix_normal_scan [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix operator inv_op (identity 6, replicate 6 0f64) t

-- ==
-- entry: test_rootfix_normal_scan_C
-- script input { complex_scan_input_C 100i64 }  
-- script input { complex_scan_input_C 1000i64 }  
-- script input { complex_scan_input_C 10000i64 }  
-- script input { complex_scan_input_C 100000i64 }  
-- script input { complex_scan_input_C 1000000i64 }  
-- script input { complex_scan_input_C 2000000i64 }  
-- script input { complex_scan_input_C 4000000i64 }  
-- script input { complex_scan_input_C 8000000i64 }  
entry test_rootfix_normal_scan_C [n] (data : [n](X_Compact, mv)) (lp : [n]i64) (rp : [n]i64) : [n](X_Compact, mv) =
    let t = T.lprp <| mkt2 lp rp data
    in T.irootfix operator_C inv_op_C (copy transform_identity, {w = [0,0,0], v_O = [0,0,0]}) t

-- ==
-- entry: test_work_efficient_rootfix_against_exscan
-- input { 256i64 }  output {true}
-- input { 4016i64 }  output {true}
-- input { 6384i64 }  output {true}
entry test_work_efficient_rootfix_against_exscan (n : i64) : bool =
    let (data, lp, rp) = complex_scan_input n
    let t = T.lprp <| mkt2 lp rp data
    let a = T.rootfix_blocked operator inv_op (identity 6, replicate 6 0f64) t 256i64
    let b = T.rootfix operator inv_op (identity 6, replicate 6 0f64) t
    in map2 (\xr yr -> 
                map2 (\x y -> f64.abs (x - y) < 1e-8f64) xr yr
                |> reduce (&&) true
            ) (map (.1) a) (map (.1) b)
        |> and



-- ==
-- entry: bench_blocked
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32} 
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry bench_blocked [n] (xs : [n]i32) : [n]i32 =
    blocked (+) 0 xs 64

-- ==
-- entry: bench_work_efficient
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32} 
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry bench_work_efficient [n] (xs : [n]i32) : [n]i32 =
    work_efficient (+) 0 xs 

-- ==
-- entry: bench_regular_scan
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32} 
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry bench_regular_scan [n] (xs : [n]i32) : [n]i32 =
    scan (+) 0 xs

