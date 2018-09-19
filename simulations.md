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
consense -f ExaBayes_topologies.run-0.024 -n 024.cons
## Printed to screen:
##Printed consensus tree in nexus format to ExaBayes_ConsensusExtendedMajorityRuleNexus.024.cons
##Printed consensus tree in newick format to ExaBayes_ConsensusExtendedMajorityRuleNewick.024.cons
consense -f ExaBayes_topologies.run-0.027 -n 027.cons
consense -f ExaBayes_topologies.run-0.036 -n 036.cons
consense -f ExaBayes_topologies.run-0.043 -n 043.cons
consense -f ExaBayes_topologies.run-0.050 -n 050.cons
consense -f ExaBayes_topologies.run-0.064 -n 064.cons
consense -f ExaBayes_topologies.run-0.071 -n 071.cons
consense -f ExaBayes_topologies.run-0.125 -n 125.cons
consense -f ExaBayes_topologies.run-0.150 -n 150.cons
consense -f ExaBayes_topologies.run-0.354 -n 354.cons
consense -f ExaBayes_topologies.run-0.404 -n 404.cons
consense -f ExaBayes_topologies.run-0.500 -n 500.cons
```

- Credible set for trees:
```shell
credibleSet -f ExaBayes_topologies.run-0.024 -n 024.cred
##printed the 50% credible set to ExaBayes_credibleSet.024.cred
credibleSet -f ExaBayes_topologies.run-0.027 -n 027.cred
credibleSet -f ExaBayes_topologies.run-0.036 -n 036.cred
credibleSet -f ExaBayes_topologies.run-0.043 -n 043.cred
credibleSet -f ExaBayes_topologies.run-0.050 -n 050.cred
credibleSet -f ExaBayes_topologies.run-0.064 -n 064.cred
credibleSet -f ExaBayes_topologies.run-0.071 -n 071.cred
credibleSet -f ExaBayes_topologies.run-0.125 -n 125.cred
credibleSet -f ExaBayes_topologies.run-0.150 -n 150.cred
credibleSet -f ExaBayes_topologies.run-0.354 -n 354.cred
credibleSet -f ExaBayes_topologies.run-0.404 -n 404.cred
credibleSet -f ExaBayes_topologies.run-0.500 -n 500.cred
```

- Parameters:
```shell
postProcParam -f ExaBayes_parameters.run-0.024 -n 024.pars
#reading file ExaBayes_parameters.run-0.024
#successfully printed statistics to file >ExaBayes_parameterStatistics.024.pars<
postProcParam -f ExaBayes_parameters.run-0.027 -n 027.pars
postProcParam -f ExaBayes_parameters.run-0.036 -n 036.pars
postProcParam -f ExaBayes_parameters.run-0.043 -n 043.pars
postProcParam -f ExaBayes_parameters.run-0.050 -n 050.pars
postProcParam -f ExaBayes_parameters.run-0.064 -n 064.pars
postProcParam -f ExaBayes_parameters.run-0.071 -n 071.pars
postProcParam -f ExaBayes_parameters.run-0.125 -n 125.pars
postProcParam -f ExaBayes_parameters.run-0.150 -n 150.pars
postProcParam -f ExaBayes_parameters.run-0.354 -n 354.pars
postProcParam -f ExaBayes_parameters.run-0.404 -n 404.pars
postProcParam -f ExaBayes_parameters.run-0.500 -n 500.pars
```
Alternatively, you could also open the parameter file with `Tracer` and visualize the distributions. If you do not have `Tracer` installed, have a look at `ExaBayes_parameterStatistics.params` (spreadsheet tools like Excel are helpful). You'll find summary statistics for each parameter. Specifically, check out the effective sampling size (ESS) value for each parameter. We will use this later when comparing to bistro ESS.

3. Use the consensus tree to simulate data with our script: `Bistro/Scripts/simulatingData-paper.r` which is the newer script (nsites=1500). We need the consensus tree from exabayes, pi and r.
```shell
(master) $ head ExaBayes_ConsensusExtendedMajorityRuleNewick.024.cons
(Cassowary:3.615832728653814e-02,Emu:4.154326867301236e-02,((BrownKiwi:1.667964241507583e-02,(LSKiwi:6.111469276818545e-03,gskiwi:5.616096315049349e-03)1.00:1.338927883291753e-02)1.00:8.038596520936106e-02,(((Dinornis:1.223356897152369e-02,(EasternMoa:6.732163590818093e-03,lbmoa:6.459366117859303e-03)1.00:8.069998891244120e-03)1.00:7.923690992359837e-02,(ECtinamou:1.878516508710195e-01,(Gtinamou:1.462000904222475e-01,Crypturellus:2.090189407003211e-01)1.00:4.876373915807834e-02)1.00:1.062879087015135e-01)1.00:3.031823931399621e-02,((LesserRhea:1.704142995190647e-02,GreatRhea:1.700795659411470e-02)1.00:1.146965660062112e-01,(Ostrich:1.214248466520340e-01,((Alligator:1.822883900249740e-01,Caiman:3.670001271956478e-01)1.00:1.750960883748059e+00,(((magpiegoose:9.347735943733071e-02,duck:1.761784440111214e-01)1.00:2.472232098666592e-02,(Chicken:1.813695928474492e-01,BrushTurkey:1.406923408018491e-01)1.00:4.088480521967741e-02)1.00:4.493295486484707e-02,((GaviaStellata:8.917191869582321e-02,LBPenguin:1.181315066517599e-01)1.00:1.347249009571150e-02,(turnstone:1.045094806127274e-01,oystercatcher:1.078915485422830e-01)1.00:1.869964773559406e-02)1.00:4.693012403498718e-02)1.00:5.499630553295962e-02)1.00:5.170227475238045e-02)0.64:1.036568160155685e-02)0.67:1.034432492772586e-02)1.00:1.088540042721619e-02)1.00:4.644965531543793e-02);
(master) $ head ExaBayes_parameterStatistics.024.pars
paramName	mean	sd	perc5	perc25	median	perc75	per95	ESS
r{0}(G<->T)	5.680462272613145e-02	4.613515885982771e-03	5.660941440021038e-02	5.368855941508072e-02	5.660941440021038e-02	5.992894121127219e-02	6.439271061477057e-02	6.443279500271909e+01
r{0}(C<->T)	3.660462779068075e-01	9.769555798416069e-03	3.657608912331709e-01	3.596325675088620e-01	3.657608912331709e-01	3.722620162793180e-01	3.827763698450830e-01	3.129634942579820e+02
r{0}(A<->G)	3.253110631174742e-01	9.183309123101020e-03	3.254830995765860e-01	3.191811007048929e-01	3.254830995765860e-01	3.314407602423361e-01	3.400639421652209e-01	3.883080897274309e+02
alpha{0}	2.335720505097140e-01	3.662125709475038e-03	2.335614639310523e-01	2.310643552857518e-01	2.335614639310523e-01	2.361270537913010e-01	2.395042423666564e-01	1.125341546754005e+03
pi{0}(G)	1.635302814698019e-01	2.933800114214484e-03	1.633902319297486e-01	1.616579548633981e-01	1.633902319297486e-01	1.655294952235767e-01	1.682106402661479e-01	5.001445658057845e+02
pi{0}(T)	2.637257253897267e-01	3.309941387192296e-03	2.637421939848149e-01	2.614534349043434e-01	2.637421939848149e-01	2.659339195542327e-01	2.694255169800950e-01	5.843523336173172e+02
r{0}(A<->C)	1.023267905961401e-01	4.480793426053360e-03	1.021612882292625e-01	9.928102458583188e-02	1.021612882292625e-01	1.053777917833203e-01	1.098375259192818e-01	4.296094716622378e+02
pi{0}(C)	2.899478977685102e-01	3.316896309610341e-03	2.900350901385625e-01	2.877633050328054e-01	2.900350901385625e-01	2.922083928054798e-01	2.953324504197168e-01	6.827055679446304e+02
r{0}(C<->G)	3.444293045562889e-02	3.386263921146780e-03	3.436077708567455e-02	3.213892131079511e-02	3.436077708567455e-02	3.663144823716875e-02	4.007994556183447e-02	8.300317659375557e+02
```

Now, we run the R script that will simulate data on the exabayes tree, rates vector r and pi vector:
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Scripts
Rscript --vanilla simulatingData-paper.r 024 1500 1633
Rscript --vanilla simulatingData-paper.r 027 1500 1952
Rscript --vanilla simulatingData-paper.r 036 1500 4007
Rscript --vanilla simulatingData-paper.r 043 1500 8300
Rscript --vanilla simulatingData-paper.r 050 1500 3213
Rscript --vanilla simulatingData-paper.r 064 1500 6827
Rscript --vanilla simulatingData-paper.r 071 1500 4296
Rscript --vanilla simulatingData-paper.r 125 1500 4480
Rscript --vanilla simulatingData-paper.r 150 1500 1023
Rscript --vanilla simulatingData-paper.r 354 1500 9928
Rscript --vanilla simulatingData-paper.r 404 1500 5843
Rscript --vanilla simulatingData-paper.r 500 1500 2637
```

This script will create a fasta file (for bistro) and a nexus file (originally for mrbayes), but we want to compare to exabayes, so we will convert the nexus file to phylip file with `scripts/nexus2phylip.nex`

```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Data/datasets
perl ../../../scripts/nexus2phylip.pl -nexus sim-024-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-027-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-036-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-043-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-050-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-064-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-071-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-125-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-150-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-354-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-404-nsites-1500.nex
perl ../../../scripts/nexus2phylip.pl -nexus sim-500-nsites-1500.nex
```

4. We run bistro on the simulated data with/without fixed tree, see `Bistro/Scripts/bistroOneRep.pl` (need to decide the sample size)

- Compile bistro
```shell
cd Documents/github/CladeCondProb/ccdprobs/Bistro/Code/bistro
cp bistro bistro0 ##backup
make ## errors, so we will use the old executable: tell bret
```

**aqui voy**: necesito modificar runningBistro.jl para correr bistro con fixed tree y sin fixed tree, hacer el script para correr un unico dataset, y correr en paralelo

5. Run exabayes on the simulated sequences (save time)
6. Compare ESS of branch lengths and Q between bistro and exabayes