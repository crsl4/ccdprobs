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
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 027
exabayes -f 027.phy -n 027 -s 1130 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 036
exabayes -f 036.phy -n 036 -s 700 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 043
exabayes -f 043.phy -n 043 -s 530 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 050
exabayes -f 050.phy -n 050 -s 4045 -m DNA
```
```shell
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
/s/std/bin/stashticket 
/s/std/bin/runauth /usr/bin/screen -S 064
exabayes -f 064.phy -n 064 -s 6509 -m DNA
```
All started on 9/15/18, 5pm. **aqui voy**

2. Summarizing exabayes results: we want to extract a consensus tree. We need Q and pi as well.
3. Use the consensus tree to simulate data with our script: `Bistro/Scripts/simulatingData.r` or `Bistro/Scripts/simulateDataStudy.r` (nsites=500,1500)
4. We run bistro on the simulated data with/without fixed tree, see `Bistro/Scripts/bistroOneRep.pl` (need to decide the sample size)
5. Run exabayes on the simulated sequences (save time)
6. Compare ESS of branch lengths and Q between bistro and exabayes