# To do
- wait for runs in mac desktop
- email mike c with error in darwin

# Simulations for paper 

0. Copy all phylip files to darwin02:
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
scp *.phy claudia@darwin02.stat.wisc.edu:/u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
```


1. We will run exabayes in 024,027,036,043,050,064. We cannot run more than 1 chain due to error (see `error-exabayes.md`), so we will run sequentially. We can only use `darwin02`. See `results.md` for a detailed description of the output files.

```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 024
exabayes -f 024.phy -n 024 -s 1234 -m DNA

exabayes: ./src/mcmc/SampleMaster.cpp:1326: void SampleMaster::writeCheckpointMaster(): Assertion `ret == 0' failed.
[darwin02:32561] *** Process received signal ***
[darwin02:32561] Signal: Aborted (6)
[darwin02:32561] Signal code:  (-6)
[darwin02:32561] [ 0] /lib/x86_64-linux-gnu/libpthread.so.0(+0x10330) [0x7f7ca4389330]
[darwin02:32561] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x37) [0x7f7ca3fe6c37]
[darwin02:32561] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x148) [0x7f7ca3fea028]
[darwin02:32561] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x2fbf6) [0x7f7ca3fdfbf6]
[darwin02:32561] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x2fca2) [0x7f7ca3fdfca2]
[darwin02:32561] [ 5] exabayes() [0x426767]
[darwin02:32561] [ 6] exabayes() [0x427a78]
[darwin02:32561] [ 7] exabayes() [0x4118f9]
[darwin02:32561] [ 8] exabayes() [0x411eba]
[darwin02:32561] [ 9] exabayes() [0x40f7ac]
[darwin02:32561] [10] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf5) [0x7f7ca3fd1f45]
[darwin02:32561] [11] exabayes() [0x41127f]
[darwin02:32561] *** End of error message ***
Aborted (core dumped)
```
We ran again in mac desktop:
```shell
[1000000,1.55s]	 -95,298.75
best state for run 0 was: -95286.13

Converged/stopped after 1000000 generations

Total walltime elapsed:	3964.59 seconds 	or 01:06:4.59 (hh:mm:ss).
Total CPU time elapsed:	3964.59 seconds 	or 01:06:4.59 (hh:mm:ss).
```

```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 027
exabayes -f 027.phy -n 027 -s 1130 -m DNA

[1000000,0.48s]  -6,514.73
best state for run 0 was: -6494.41

Converged/stopped after 1000000 generations

Total walltime elapsed: 1113.00 seconds         or 00:18:33.00 (hh:mm:ss).
Total CPU time elapsed: 1113.00 seconds         or 00:18:33.00 (hh:mm:ss).
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 036
exabayes -f 036.phy -n 036 -s 700 -m DNA

[1000000,0.57s]  -27,100.75
best state for run 0 was: -27071.39

Converged/stopped after 1000000 generations

Total walltime elapsed: 1370.86 seconds         or 00:22:50.86 (hh:mm:ss).
Total CPU time elapsed: 1370.86 seconds         or 00:22:50.86 (hh:mm:ss).
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 043
exabayes -f 043.phy -n 043 -s 530 -m DNA

exabayes: ./src/mcmc/SampleMaster.cpp:1326: void SampleMaster::writeCheckpointMaster(): Assertion `ret == 0' failed.
[darwin02:32868] *** Process received signal ***
[darwin02:32868] Signal: Aborted (6)
[darwin02:32868] Signal code:  (-6)
[darwin02:32868] [ 0] /lib/x86_64-linux-gnu/libpthread.so.0(+0x10330) [0x7fb6b5f59330]
[darwin02:32868] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x37) [0x7fb6b5bb6c37]
[darwin02:32868] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x148) [0x7fb6b5bba028]
[darwin02:32868] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x2fbf6) [0x7fb6b5bafbf6]
[darwin02:32868] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x2fca2) [0x7fb6b5bafca2]
[darwin02:32868] [ 5] exabayes() [0x426767]
[darwin02:32868] [ 6] exabayes() [0x427a78]
[darwin02:32868] [ 7] exabayes() [0x4118f9]
[darwin02:32868] [ 8] exabayes() [0x411eba]
[darwin02:32868] [ 9] exabayes() [0x40f7ac]
[darwin02:32868] [10] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf5) [0x7fb6b5ba1f45]
[darwin02:32868] [11] exabayes() [0x41127f]
[darwin02:32868] *** End of error message ***
Aborted (core dumped)
```
We ran again in mac desktop:
```shell
[1000000,0.76s]	 -17,466.61
best state for run 0 was: -17447.34

Converged/stopped after 1000000 generations

Total walltime elapsed:	1654.58 seconds 	or 00:27:34.58 (hh:mm:ss).
Total CPU time elapsed:	1654.58 seconds 	or 00:27:34.58 (hh:mm:ss).
```

```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 050
exabayes -f 050.phy -n 050 -s 4045 -m DNA

[1000000,0.97s]  -5,976.87
best state for run 0 was: -5939.79

Converged/stopped after 1000000 generations

Total walltime elapsed: 1130.02 seconds         or 00:18:50.02 (hh:mm:ss).
Total CPU time elapsed: 1130.02 seconds         or 00:18:50.02 (hh:mm:ss).
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 064
exabayes -f 064.phy -n 064 -s 6509 -m DNA

[1000000,1.02s]  -7,093.87
best state for run 0 was: -7061.68

Converged/stopped after 1000000 generations

Total walltime elapsed: 1044.05 seconds         or 00:17:24.05 (hh:mm:ss).
Total CPU time elapsed: 1044.05 seconds         or 00:17:24.05 (hh:mm:ss).
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 071
exabayes -f 071.phy -n 071 -s 453 -m DNA

Error while opening/reading file tmp.201583265

Error: parsing file tmp.201583265 failed. 

Please double-check, whether each site is assigned a partition. 
Side note: the notation for choosing every n-th (e.g., triplet) 
exabayes: ./src/contrib/AlignmentPLL.cpp:208: void AlignmentPLL::initPartitions(std::string): Assertion `0' failed.
character has changed from '\3' to '/3'.[darwin02:37024] *** Process received signal ***
[darwin02:37024] Signal: Aborted (6)
[darwin02:37024] Signal code:  (-6)
[darwin02:37024] [ 0] /lib/x86_64-linux-gnu/libpthread.so.0(+0x10330) [0x7fa208ecd330]
[darwin02:37024] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x37) [0x7fa208b2ac37]
[darwin02:37024] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x148) [0x7fa208b2e028]
[darwin02:37024] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x2fbf6) [0x7fa208b23bf6]
[darwin02:37024] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x2fca2) [0x7fa208b23ca2]
[darwin02:37024] [ 5] exabayes() [0x452b07]
[darwin02:37024] [ 6] exabayes() [0x4546ff]
[darwin02:37024] [ 7] exabayes() [0x4285d4]
[darwin02:37024] [ 8] exabayes() [0x429953]
[darwin02:37024] [ 9] exabayes() [0x4118c2]
[darwin02:37024] [10] exabayes() [0x411eba]
[darwin02:37024] [11] exabayes() [0x40f7ac]
[darwin02:37024] [12] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf5) [0x7fa208b15f45]
[darwin02:37024] [13] exabayes() [0x41127f]
[darwin02:37024] *** End of error message ***
Aborted (core dumped)
```

We will run in Mac desktop in the meantime:
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 071.phy -n 071 -s 453 -m DNA

[1000000,0.85s]	 -5,368.99
best state for run 0 was: -5330.01

Converged/stopped after 1000000 generations

Total walltime elapsed:	1278.26 seconds 	or 00:21:18.26 (hh:mm:ss).
Total CPU time elapsed:	1278.26 seconds 	or 00:21:18.26 (hh:mm:ss).
```

```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 125.phy -n 125 -s 6619 -m DNA

[1000000,10.13s]	 -824,132.63
best state for run 0 was: -824087.65

Converged/stopped after 1000000 generations

Total walltime elapsed:	26491.67 seconds 	or 07:21:31.67 (hh:mm:ss).
Total CPU time elapsed:	26491.67 seconds 	or 07:21:31.67 (hh:mm:ss).
```
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 140.phy -n 140 -s 1343 -m DNA

initialized diagnostics file ExaBayes_diagnostics.140
initialized file ExaBayes_topologies.run-0.140
initialized file ExaBayes_parameters.run-0.140
[iMac:52067] *** Process received signal ***
[iMac:52067] Signal: Segmentation fault: 11 (11)
[iMac:52067] Signal code: Address not mapped (1)
[iMac:52067] Failing at address: 0x7ffeedfca600
[iMac:52067] [ 0] 0   libsystem_platform.dylib            0x00007fff54a85f5a _sigtramp + 26
[iMac:52067] [ 1] 0   ???                                 0x0000000000000000 0x0 + 0
[iMac:52067] [ 2] 0   exabayes                            0x0000000101d1699d pllNewviewIterative + 6316
[iMac:52067] *** End of error message ***
Segmentation fault: 11
```
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 150.phy -n 150 -s 220 -m DNA

[1000000,1.77s]	 -39,770.35
best state for run 0 was: -39716.42

Converged/stopped after 1000000 generations

Total walltime elapsed:	3734.78 seconds 	or 01:02:14.78 (hh:mm:ss).
Total CPU time elapsed:	3734.78 seconds 	or 01:02:14.78 (hh:mm:ss).
```
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 354.phy -n 354 -s 1819 -m DNA

[1000000,1.46s]	 -7,264.87
best state for run 0 was: -7197.48

Converged/stopped after 1000000 generations

Total walltime elapsed:	4281.95 seconds 	or 01:11:21.95 (hh:mm:ss).
Total CPU time elapsed:	4281.95 seconds 	or 01:11:21.95 (hh:mm:ss).
```
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 404.phy -n 404 -s 1113 -m DNA

[1000000,9.95s]	 -156,595.63
best state for run 0 was: -156524.69

Converged/stopped after 1000000 generations

Total walltime elapsed:	22467.89 seconds 	or 06:14:27.89 (hh:mm:ss).
Total CPU time elapsed:	22467.89 seconds 	or 06:14:27.89 (hh:mm:ss).
```
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
exabayes -f 500.phy -n 500 -s 7704 -m DNA

[1000000,3.44s]	 -86,336.12
best state for run 0 was: -86240.38

Converged/stopped after 1000000 generations

Total walltime elapsed:	8525.69 seconds 	or 02:22:5.69 (hh:mm:ss).
Total CPU time elapsed:	8525.69 seconds 	or 02:22:5.69 (hh:mm:ss).
```

024-064 started on 9/15/18, 5pm, all finished within 20 minutes.
071-500 started on 9/16/18, 5pm, on desktop, all finished within 2 hours.

First, we will copy all results from darwin into repo:
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
scp claudia@darwin02.stat.wisc.edu:/u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes/ExaBayes_*.027 .
scp claudia@darwin02.stat.wisc.edu:/u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes/ExaBayes_*.036 .
scp claudia@darwin02.stat.wisc.edu:/u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes/ExaBayes_*.050 .
scp claudia@darwin02.stat.wisc.edu:/u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes/ExaBayes_*.064 .
```

2. Summarizing exabayes results: we want to extract a consensus tree. We need Q and pi as well.

- Consensus tree:
```shell
consense -f ExaBayes_topologies.024.0 -n 024.cons
consense -f ExaBayes_topologies.027.0 -n 027.cons
consense -f ExaBayes_topologies.036.0 -n 036.cons
consense -f ExaBayes_topologies.043.0 -n 043.cons
consense -f ExaBayes_topologies.050.0 -n 050.cons
consense -f ExaBayes_topologies.064.0 -n 064.cons
consense -f ExaBayes_topologies.071.0 -n 071.cons
consense -f ExaBayes_topologies.125.0 -n 125.cons
consense -f ExaBayes_topologies.150.0 -n 150.cons
consense -f ExaBayes_topologies.354.0 -n 354.cons
consense -f ExaBayes_topologies.404.0 -n 404.cons
consense -f ExaBayes_topologies.500.0 -n 500.cons
```

3. Use the consensus tree to simulate data with our script: `Bistro/Scripts/simulatingData.r` or `Bistro/Scripts/simulateDataStudy.r` (nsites=500,1500)
4. We run bistro on the simulated data with/without fixed tree, see `Bistro/Scripts/bistroOneRep.pl` (need to decide the sample size)
5. Run exabayes on the simulated sequences (save time)
6. Compare ESS of branch lengths and Q between bistro and exabayes