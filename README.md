# Data Parallel Rigid Body Dynamics
This repo contains data parallel code for rigid body dynamics. A key datastructure used to achieve this is the vector tree defined by Blelloch in [Ble90](https://www.cs.cmu.edu/~guyb/papers/Ble90.pdf) p. 84-91.


### Matlab On Server

Check out this article to see how to set up matlab on a server and run scripts without a GUI:
<https://www.funwithlinux.net/blog/run-matlab-in-linux-without-graphical-environment/>

The matlab code that is tested against is made by Roy Featherstone and can be found here: <https://royfeatherstone.org/spatial/v2/index.html>

### Benchmarks

The benchmarks for RNEA inside the bench folder can be run with:

```
futhark bench  bench_rnea.fut
```
If you want to benchmark using the work efficient scan when doing rootfix and leaffix operation you need to replace the vtree.fut file from the vtree futhark package with the "vtree_with_work_efficient_scan.fut". E.g. with:

```bash
cp lib/github.com/diku-dk/vtree/vtree.fut lib/github.com/diku-dk/vtree/vtree_backup.fut 
rm lib/github.com/diku-dk/vtree/vtree.fut 
cp vtree_with_work_efficient_scan.fut lib/github.com/diku-dk/vtree/vtree.fut 
```


