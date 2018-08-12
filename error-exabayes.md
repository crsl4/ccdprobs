# Using exabayes example data
Following section 4.2 in [manual](https://sco.h-its.org/exelixis/web/software/exabayes/manual/index.html#comp-src), we get an error (not very informative):
```shell
[claudia@darwin02] (40)$ yggdrasil -f aln.phy -q aln.part  -n myRun -s $RANDOM -c config.nex

Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!

[claudia@darwin02] (41)$ mpirun -np 16 exabayes -f aln.phy -q aln.part  -n myRun -s $RANDOM -c config.nex -R 4

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD 
with errorcode -1.

NOTE: invoking MPI_ABORT causes Open MPI to kill all MPI processes.
You may or may not see output from other processes, depending on
exactly when Open MPI kills them.
--------------------------------------------------------------------------
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
Error: tried to create intermediate file ExaBayes_binaryAlignment.myRun, but did not succeed!
--------------------------------------------------------------------------
mpirun has exited due to process rank 2 with PID 32147 on
node darwin02 exiting improperly. There are two reasons this could occur:

1. this process did not call "init" before exiting, but others in
the job did. This can cause a job to hang indefinitely while it waits
for all processes to call "init". By rule, if one process calls "init",
then ALL processes must call "init" prior to termination.

2. this process called "init", but exited without calling "finalize".
By rule, all processes that call "init" MUST call "finalize" prior to
exiting or it will be considered an "abnormal termination"

This may have caused other processes in the application to be
terminated by signals sent by mpirun (as reported here).
--------------------------------------------------------------------------
[darwin02:32144] 15 more processes have sent help message help-mpi-api.txt / mpi-abort
[darwin02:32144] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
```

# Using our own data
But this could be a problem with the example data.
So, let's try to run with our data (something that worked already with `yggdrasil`):
```shell
ssh claudia@darwin02.stat.wisc.edu
cd /u/c/l/claudia/Documents/phylo/projects/CladeDist/exabayes
mpirun -np 2 exabayes -f 024.phy -n 1 -s 1234 -m DNA -R 2

This is the multi-threaded MPI hybrid variant of ExaBayes (version 1.5),
a tool for Bayesian MCMC sampling of phylogenetic trees, build with the
Phylogenetic Likelihood Library (version 1.0.0, September 2013).

This software has been released on 2016-09-24 13:24:08
(git commit id:ef87cd16821f979075dcc20c870298d87406f146)

	by Andre J. Aberer, Kassian Kobert and Alexandros Stamatakis

Please send any bug reports, feature requests and inquiries to exabayes-at-googlegroups-dot-com

The program was called as follows: 
exabayes -f 024.phy -n 1 -s 1234 -m DNA -R 2 

================================================================
You provided an alignment file in phylip format. Trying to parse it...

	ERROR: you want to run  2 independent runs in parallel, while there are only 1 to be run in total.


--------------------------------------------------------------------------
mpirun noticed that the job aborted, but has no info as to the process
that caused that situation.
--------------------------------------------------------------------------
```

New test with option `-C` instead of `-R`:
```shell
mpirun -np 2 exabayes -f 024.phy -n 1 -s 1234 -m DNA -C 2


This is the multi-threaded MPI hybrid variant of ExaBayes (version 1.5),
a tool for Bayesian MCMC sampling of phylogenetic trees, build with the
Phylogenetic Likelihood Library (version 1.0.0, September 2013).

This software has been released on 2016-09-24 13:24:08
(git commit id:ef87cd16821f979075dcc20c870298d87406f146)

	by Andre J. Aberer, Kassian Kobert and Alexandros Stamatakis

Please send any bug reports, feature requests and inquiries to exabayes-at-googlegroups-dot-com

The program was called as follows: 
exabayes -f 024.phy -n 1 -s 1234 -m DNA -C 2 

================================================================
You provided an alignment file in phylip format. Trying to parse it...

	ERROR: you want to run 2 coupled chains in parallel, while there are only 1 to be run in total.


--------------------------------------------------------------------------
mpirun noticed that the job aborted, but has no info as to the process
that caused that situation.
--------------------------------------------------------------------------
```

Finally, not using multiple chains (which is identical to `yggdrasil`):
```shell
exabayes -f 024.phy -n 1 -s 1234 -m DNA
## works!
```
**Note** that I get the same errors when I run in my Mac.