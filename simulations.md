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
**aqui voy** this is weird, because I get error with all the files, even files that ran before (like 064). I need to email Mike.

```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 125
exabayes -f 125.phy -n 125 -s 6619 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 140
exabayes -f 140.phy -n 140 -s 1343 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 150
exabayes -f 150.phy -n 150 -s 220 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 354
exabayes -f 354.phy -n 354 -s 1819 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 404
exabayes -f 404.phy -n 404 -s 1113 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 500
exabayes -f 500.phy -n 500 -s 7704 -m DNA
```

024-064 started on 9/15/18, 5pm, finished within 20 minutes.
071-500 started on 9/16/18, 5pm **aqui voy**

2. Summarizing exabayes results: we want to extract a consensus tree. We need Q and pi as well.
3. Use the consensus tree to simulate data with our script: `Bistro/Scripts/simulatingData.r` or `Bistro/Scripts/simulateDataStudy.r` (nsites=500,1500)
4. We run bistro on the simulated data with/without fixed tree, see `Bistro/Scripts/bistroOneRep.pl` (need to decide the sample size)
5. Run exabayes on the simulated sequences (save time)
6. Compare ESS of branch lengths and Q between bistro and exabayes