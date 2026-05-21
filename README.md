# Data Parallel Rigid Body Dynamics
This repo contains data parallel code for rigid body dynamics. A key datastructure used to achieve this is the vector tree defined by Blelloch in [Ble90](https://www.cs.cmu.edu/~guyb/papers/Ble90.pdf) p. 84-91.

### Benchmarks
Before running benchmarks you need to download some futhark libraries with:
```
futhark pkg sync
```

The benchmarks for RNEA inside the bench folder can be run with the `futhark bench` tool.
If you e.g. want to run the rnea benchmark with the CUDA backend you can run:
```
futhark bench --backend=cuda bench_rnea2.fut
```

The benchmarks in the report were run with sbatch on the Hendrix cluster. 
Before the sbatchscript can be run you need to install Futhark and CUDA. 
This can be done e.g. with the `module load` command. Remember to also insert the path to the CUDA version you use in your .bashrc

When Futhark and CUDA is installed the RNEA and CRBA benchmarks can be run with the CUDA backend by:
```
cd bench
sbatch sbatchscript.sh
```
This also runs some sequential implementations with the c backend. The results are then found in a file specified by the sbatchscript. If you want to replicate the results from the report exactly you should make sure that the inputs specified in the benchmark files are the same.

If you want to run the benchmarks seperately you can also take inspiration from the `sbatchscript.sh` file to see how you run the relevant benchmark files.


#### Scan Variations
The implementations uses a custom vtree implementation found in `vtree_with_work_efficient_scan.fut`.
The only difference from the official vtree library is that some variations on the scan algorithm are provided.
These were shown to be faster when compared to the normal scan on complex input operators.

The scan variations can be benchmarked inside the `bench/scan_benches` folder. If you want to run the replicate the results from the report you can run the `sbatchscript.sh` file inside that folder:
```
cd bench/scan_benches
sbatch sbatchscript.sh
```

### Tests
The tests are found inside the test folder can be run with the `futhark test` tool. If you e.g. want to run the CRBA tests you can do the following:
```
cd test
futhark test crba_test_against_RF_matlab.fut 
```


### Matlab


The matlab code that is tested against is made by Roy Featherstone and can be found here: <https://royfeatherstone.org/spatial/v2/index.html>


If you want to run Featherstone's matlab implmentation you need to import his library. 
You can then use the scripts `RNEA_bench.m` or `CRBA_bench.m` in the `matlab_conversion_scripts` folder to run some benchmarks of Fethearstone's implementation.  
The benchmark scripts takes an array of input sizes `Ns` and a single branchinfactor `bf`. The scripts also write the inputs to files depending on the `filename_ending` input variable.
