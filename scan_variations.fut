import "matrix_ops"
import "spatial_ops"
import "treeModel"
import "lib/github.com/diku-dk/vtree/vtree"
-- The scan implementations are my DPP implementation where I modified them to take
-- an operator as input
def ilog2 (x: i64 ) = 63 - i64.i32 (i64.clz x)

def next_power_of_two (n : i64) : i64 =
  loop acc = 2 while acc < n do 
    acc * 2

def hillis_steele [n] 'a ( xs : [n]a) (op : a -> a -> a) : [n]a =
    let m = ilog2 n
    in loop xs = copy xs for d in (iota (m )) do
        let offset =  (2 ** (d))
        let indx_to_update = iota ((n) - offset) |>  map (\x -> x + offset)
        let update_vals = map (\i -> op xs[i] xs[i - offset]) indx_to_update
        in scatter (xs) indx_to_update update_vals

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

def work_efficient_sc [n] 'a (op : a -> a -> a) (ne : a)  (xs : [n]a)  : [n]a =
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
    let final_elem = copy upswept[k-1]
    let upswept' = upswept with [k-1] = ne
    let downswept =
        loop xs' = copy (upswept') for d in m..m-1...1 do
            let offset = 2 ** (d - 1)
            let indx_to_update = iota (k/(offset)) 
                                    |> map (\x -> (x + 1) * offset - 1)
            let update_vals = 
                map2 
                    (\itu i -> if i % 2 == 1 then op xs'[itu] xs'[itu - (offset)]
                                             else xs'[itu + offset])
                    indx_to_update (indices indx_to_update)
            in (scatter xs' indx_to_update update_vals)
    let org_arr = rotate 1 downswept 
                    |> take n 
    in org_arr with [n-1] = final_elem


-- taken from https://futhark-lang.org/blog/2019-04-10-what-is-the-minimal-basis-for-futhark.html
def scan' [n] 'a (op: a -> a -> a) (_ne: a) (as: [n]a): [n]a =
  let iters = i64.f32 (f32.ceil (f32.log2 (f32.i64 n)))
  in loop as for i < iters do
       map (\j -> if j < 2**i
                  then as[j]
                  else as[j] `op` as[j-2**i])
           (iota n)


def rootfix_hillis_steele 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = scan' op ne R
    in map (\i -> S[i]) lp

def rootfix_work_efficient 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = work_efficient op ne R 
    in map (\i -> S[i]) lp

def rootfix_work_efficient_sc 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = work_efficient_sc op ne R
    in map (\i -> S[i]) lp

def rootfix2 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = scan op ne R
    in map (\i -> S[i]) lp

def rootfix_exscan 'a [n] (op: a -> a -> a) (inv: a -> a) (ne: a) (lp : [n]i64) (rp : [n]i64) (data : [n]a) : [n]a =
    let I = replicate (2 * n) ne
    let L = scatter I lp data
    let R = scatter L rp (map inv data)
    let S = exscan op ne R
    in map (\i -> S[i]) lp

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


-- ==
-- entry: test_work_efficient2
-- compiled random input { [256]i32 }  output {true}
-- compiled random input { [500]i32}   output {true}
-- compiled random input { [65536]i32}  output {true}
-- compiled random input { [104999]i32}  output {true}
entry test_work_efficient2 [n] (xs : [n]i32) : bool =
    let a = work_efficient (+) 0 xs 
    let b = exscan (+) 0 xs
    in map2 (\x y ->  (x - y) == 0) a b 
       |> and


-- ==
-- entry: test_work_efficient_rootfix_against_exscan
-- input { 256i64 }  output {true}
-- input { 4016i64 }  output {true}
-- input { 6384i64 }  output {true}
-- input { 131072i64 }  output {true}
entry test_work_efficient_rootfix_against_exscan (n : i64) : bool =
    let (data, lp, rp) = complex_scan_input n
    let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      let inv_cia = XBtoA_from_XAtoB_M ci.0
      let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
      in (inv_cia, inv_cib)
    let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
    let a = rootfix_work_efficient operator inv_op (identity 6, replicate 6 0f64) lp rp data
    let b = rootfix_exscan operator inv_op (identity 6, replicate 6 0f64) lp rp data
    in map2 (\xr yr -> 
                map2 (\x y -> f64.abs (x - y) < 1e-8f64) xr yr
                |> reduce (&&) true
            ) (map (.1) a) (map (.1) b)
        |> and

-- ==
-- entry: test_hillis_steele_rootfix
-- script input { complex_scan_input 256i64 }  
-- script input { complex_scan_input 4096i64 }  
-- script input { complex_scan_input 16384i64 }  
-- script input { complex_scan_input 131072i64 }  
entry test_hillis_steele_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      let inv_cia = XBtoA_from_XAtoB_M ci.0
      let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
      in (inv_cia, inv_cib)
    let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
    in rootfix_hillis_steele operator inv_op (identity 6, replicate 6 0f64) lp rp data

-- ==
-- entry: test_work_efficient_rootfix
-- script input { complex_scan_input 256i64 }  
-- script input { complex_scan_input 4096i64 }  
-- script input { complex_scan_input 16384i64 }  
-- script input { complex_scan_input 131072i64 }  
entry test_work_efficient_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      let inv_cia = XBtoA_from_XAtoB_M ci.0
      let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
      in (inv_cia, inv_cib)
    let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
    in rootfix_work_efficient operator inv_op (identity 6, replicate 6 0f64) lp rp data

-- ==
-- entry: test_rootfix
-- script input { complex_scan_input 256i64 }  
-- script input { complex_scan_input 4096i64 }  
-- script input { complex_scan_input 16384i64 }  
-- script input { complex_scan_input 131072i64 }  
entry test_rootfix [n] (data : [n]([6][6]f64, [6]f64)) (lp : [n]i64) (rp : [n]i64) : [n]([6][6]f64, [6]f64) =
    let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      let inv_cia = XBtoA_from_XAtoB_M ci.0
      let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
      in (inv_cia, inv_cib)
    let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
      (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
    in rootfix2 operator inv_op (identity 6, replicate 6 0f64) lp rp data

-- ==
-- entry: test_hillis_steele
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32}  auto output
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry test_hillis_steele [n] (xs : [n]i32) : [n]i32 =
    hillis_steele xs (+)

-- ==
-- entry: test_work_efficient
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32} 
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry test_work_efficient [n] (xs : [n]i32) : [n]i32 =
    work_efficient (+) 0 xs 

-- ==
-- entry: test_regular_scan
-- compiled random input { [256]i32 }  auto output
-- compiled random input { [4096]i32} 
-- compiled random input { [65536]i32}  auto output
-- compiled random input { [1048576]i32}  auto output
-- compiled random input { [4194304]i32}  auto output
entry test_regular_scan [n] (xs : [n]i32) : [n]i32 =
    scan (+) 0 xs


-- def main =
    -- let n = 4
    -- let (data, lp, rp) = complex_scan_input n
    -- let inv_op (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    --   let inv_cia = XBtoA_from_XAtoB_M ci.0
    --   let inv_cib = scal_mul_vec_f64 (-1) (mat_mul_vec_f64 inv_cia ci.1)
    --   in (inv_cia, inv_cib)
    -- let operator (si : ([6][6]f64, [6]f64)) (ci : ([6][6]f64, [6]f64)) : ([6][6]f64, [6]f64) =
    --   (ci.0 `matmul_f64` si.0,    (ci.0 `mat_mul_vec_f64` si.1) `vecadd_f64` ci.1)
    -- let a = rootfix_work_efficient operator inv_op (identity 6, replicate 6 0f64) lp rp data
    -- let b = rootfix_exscan         operator inv_op (identity 6, replicate 6 0f64) lp rp data
    -- in (ammm, bmmm)

    -- let xs = [3i32, 1,7, 0, 4]
    -- let a = trace <| work_efficient (+) 0 xs 
    -- let b = trace <| exscan (+) 0 xs
    -- in map2 (\x y ->  (x - y) == 0) a b 
    --    |> and
