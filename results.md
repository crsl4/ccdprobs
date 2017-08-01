# Results

## Simulated whales

- When we fix the topology to the true tree (or mean tree), we get an `ESS 8%`
- However, when we do not fix the topology, we get an `ESS 1.7%`
- We want to study now the topology sampling

```r
source("../../Scripts/summary.R")
compareBistroMB("bistro4-whales",truetree="(1:0.0231136,2:0.0372603,((3:0.094114,(((4:0.025742,5:0.022279):0.0249355,14:0.0685793):0.00162364,((6:0.027438,7:0.0173734):0.0210568,((((8:0.0170381,9:0.0189656):0.0173591,11:0.0375513):0.00117174,(10:0.0324844,12:0.0337327):0.00110449):0.000772143,13:0.0327726):0.00476205):0.00772679):0.014146):0.00534651,(15:0.101767,((((16:0.0503164,17:0.0470883):0.0368309,((18:0.0110323,19:0.00505243):0.02207,(20:0.0116344,21:0.0071754):0.0214576):0.0500723):0.0548426,((24:0.0799893,((25:0.0553201,26:0.0455578):0.0283923,((27:0.0766089,28:0.0689195):0.00240495,(29:0.0489559,30:0.0222969):0.0290598):0.00387214):0.00247192):0.0113158,31:0.12489):0.0064916):0.00406787,(22:0.0938027,23:0.0901243):0.00790296):0.00718762):0.0331146):0.064614);",trueq=c(0.231333,0.225334,0.173454,0.369879,0.151508,0.198957,0.0548603,0.0118846,0.578106,0.00468519))
[1] "Most frequently sampled tree is NOT equal to the true tree: (1,2,((3,(((4,5),14),((6,7),(((8,9),((10,12),13)),11)))),(15,((((16,17),((18,19),(20,21))),(((24,(25,26)),((27,28),(29,30))),31)),(22,23))))); 0.02"
[1] "true tree sampled: 0 times"
```
The true tree is:
`(1,2,((3,(((4,5),14),((6,7),((((8,9),11),(10,12)),13)))),(15,((((16,17),((18,19),(20,21))),((24,((25,26),((27,28),(29,30)))),31)),(22,23)))));`
The mean tree is:
`(1,2,((3,(((4,5),14),((6,7),(((8,9),((10,12),13)),11)))),(15,((((16,17),((18,19),(20,21))),(((24,(25,26)),((27,28),(29,30))),31)),(22,23)))));`
MrBayes MAP tree (PP=0.108054), not the same as true tree:
`(2,((((23,22),((31,(((30,29),(28,27)),((26,25),24))),(((21,20),(19,18)),(17,16)))),15),(((((13,(12,10)),(11,(9,8))),(7,6)),(14,(5,4))),3)),1);`

The only difference between the true tree and the mean tree is the following clade, basically just the placement of 11 (or 13):
- True tree: `((((8,9),11),(10,12)),13)`
- Mean tree: `(((8,9),((10,12),13)),11)`
- MrBayes tree: `((13,(12,10)),(11,(9,8)))`

The true tree was sampled 9 times in the bootstrap trees. The distance between the true tree and the mean tree is (`bistro4-whales.distances`):
```
12 0.0500342
310 0.0360318
312 0.0342358
435 0.0393727
455 0.0393564
545 0.0440405
793 0.0361467
927 0.0376625
946 0.0412214
```
```r
> tab = read.table("bistro4-whales.distances",header=FALSE)
> str(tab)
'data.frame':	1000 obs. of  2 variables:
 $ V1: Factor w/ 750 levels "(1,(2,((3,((((15,22),(((24,(25,26)),(27,(28,(29,30)))),31)),23),((16,17),((18,19),(20,21))))),(((4,5),14),((6,7),((8,9),(((10,1"| __truncated__,..: 733 505 543 630 435 533 562 256 668 561 ...
 $ V2: num  0.0476 0.0436 0.0386 0.0425 0.0455 ...
 > summary(tab$V2)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.03019 0.03851 0.04190 0.04208 0.04537 0.05942
> ind=c(12,310,312,435,455,545,793,927,946)
> summary(tab$V2[ind])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.03424 0.03649 0.03936 0.03984 0.04122 0.05003
```
So, the true tree has on average slightly smaller distance than all the bootstrap trees. Why is it never sampled in the importance sampling?

In `bistro4-whales.topCounts`, we have a table with the bootstrap trees, with colums: count parsimonyWt parsimonyScore parsimonyDiff distanceWt
```
(1,2,((3,(((4,5),14),((6,7),((((8,9),11),(10,12)),13)))),(15,((((16,17),((18,19),(20,21))),((24,((25,26),((27,28),(29,30)))),31)),(22,23)))));     1     0.1353  2621   -4     0.0002
```
The normalized distance weight for the true tree is `0.0005681818`, so no surprise that we do not sample it `1000*0.0005681818=0.5681818`:
```r
> tab2 = read.table("bistro4-whales.topCounts",header=TRUE)
> str(tab2)
'data.frame':	750 obs. of  6 variables:
 $ tree          : Factor w/ 750 levels "(1,2,((3,((((15,22),(((24,(25,26)),(27,(28,(29,30)))),31)),23),((16,17),((18,19),(20,21))))),(((4,5),14),((6,7),((8,9),(((10,12"| __truncated__,..: 1 2 3 4 5 6 7 8 9 10 ...
 $ count         : int  1 1 1 1 1 1 1 1 1 1 ...
 $ parsimonyWt   : num  0e+00 0e+00 1e-04 0e+00 0e+00 0e+00 1e-04 0e+00 0e+00 0e+00 ...
 $ parsimonyScore: int  2643 2640 2636 2640 2637 2637 2636 2638 2640 2641 ...
 $ parsimonyDiff : int  -26 -23 -19 -23 -20 -20 -19 -21 -23 -24 ...
 $ distanceWt    : num  5e-04 1e-04 1e-04 0e+00 3e-04 2e-04 1e-04 1e-04 2e-04 0e+00 ...
> w=tab2$distanceWt/sum(tab2$distanceWt)
> w[126]
[1] 0.0005681818
```
The tree with biggest weight is not the mean tree:
```r
> summary(w)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
0.0000000 0.0002841 0.0008523 0.0013330 0.0014200 0.0321000
> ind2=which(w>0.031)
> ind2
[1] 562
> tab2$tree[ind2]
[1] (1,2,((3,(15,((((16,17),((18,19),(20,21))),(((24,(25,26)),((27,28),(29,30))),31)),(22,23)))),(((4,5),14),((6,7),(((8,9),((10,12),13)),11)))));
750 Levels: (1,2,((3,((((15,22),(((24,(25,26)),(27,(28,(29,30)))),31)),23),((16,17),((18,19),(20,21))))),(((4,5),14),((6,7),((8,9),(((10,12),13),11)))))); ...
> tre3=read.tree(text="(1,2,((3,(15,((((16,17),((18,19),(20,21))),(((24,(25,26)),((27,28),(29,30))),31)),(22,23)))),(((4,5),14),((6,7),(((8,9),((10,12),13)),11)))));")
> mtre=read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),((10,12),13)),11)))),(15,((((16,17),((18,19),(20,21))),(((24,(25,26)),((27,28),(29,30))),31)),(22,23)))));")
> layout(matrix(c(1,2),nrow=1))
> plot(tre3,no.margin=TRUE)
> plot(mtre,no.margin=TRUE)
```

- The true tree has the clades
  - `((8,9),11)`
  - `(((8,9),11),(10,12))`
- The mean tree has the clades
  - `((10,12),13)`
  - `((8,9),((10,12),13))`

| Clade | dist | parsimony | count | mrbayes |
|:--|:--:|:--:|:--:|
| {8-9,11} | 0.1617777 | 0.34876743 | 0.15700000 | 0.3817950 |
| {8-12} | 0.02284903 | 0.02213827 | 0.03100000 |  |
| {10,12-13} | 0.62394715 | 0.75848921 | 0.60700000 | 0.7685898 |
| {8-10,12-13} | 0.59945013 | 0.48872218 | 0.58600000 | 0.3156172 |

I think that the problem with this dataset is the very short branches around 10,11,12,13. Maybe this is a not a problem of Bistro.

## Cats and dogs
- Cats and dogs with fixed mean tree: `ESS 16%`
- not fixed tree: `ESS 2.4%` (without choosing the root wisely for the long branch)

Mean tree:
`(1,2,((3,((4,5),6)),(((((7,8),9),10),11),12)));`

Bistro tree:
```
(1,2,((3,((4,5),6)),(((((7,8),9),10),12),11))); 0.4286
(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6))); 0.3877
(1,(2,((((7,8),9),10),(11,12))),(3,((4,5),6))); 0.1191
```
MrBayes tree:
```
(((11,(12,(10,(9,(8,7))))),2),((6,(5,4)),3),1); 0.786589
(2,((11,(12,(10,(9,(8,7))))),((6,(5,4)),3)),1); 0.140961
```
The mean (bistro) tree is the second mb tree; the second bistro tree is the first mb tree. But the 1st tree in bistro matches the tree

### Simulated cats and dogs
- sim cats dogs 1500: `ESS 2.76%`
- sim cats dogs fixed tree: `ESS 16.9%`

```r
source("../../Scripts/summary.R")
compareBistroMB("bistro3-cats",mb="sim-cats-dogs-nsites1500",truetree="(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.01190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.00840558068745050208,(gray_wolf:0.00206050882985083861,dog:0.00185256446369396789):0.03205946058703370433):0.02609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;", trueq=c(0.242929500,0.235609700,0.206743300,0.314717500,0.047201519,0.293928352,0.068691153,0.012274385,0.569989255,0.007915336))
[1] "Most frequently sampled tree is NOT equal to the true tree: (1,((((2,3),5),4),(6,(((7,(8,9)),10),11))),12); 0.2"
[1] "true tree sampled: 174 times"
```
True tree:
`(1,((((2,3),4),5),(6,(((7,(8,9)),10),11))),12);`
Mean tree:
`(1,((((2,3),4),5),(6,(((7,(8,9)),10),11))),12);`

No problems here. All plots look fine (`bistro3-cats*.pdf`).

## Simulated 16 taxa
- not on fixed tree: `ESS 6.8%`
- on fixed mean tree: `ESS 13%`
```r
source("../../Scripts/summary.R")
compareBistroMB("bistro3",mb="sim-16-nsites1500",truetree="(A:0.0025,B:0.0025,((C:0.0025,D:0.0025):0.0025,(((E:0.0125,F:0.0125):0.0125,(G:0.0125,H:0.0125):0.0125):0.0625,(((((((I:0.0125,J:0.0125):0.0125,K:0.025):0.0125,L:0.0375):0.0125,M:0.05):0.0125,N:0.0625):0.0125,O:0.075):0.0125,P:0.0875):0.0125):0.0325):0.0025);",trueq=c(0.3,0.3,0.1,0.3,0.1,0.3,0.1,0.1,0.3,0.1))
[1] "Most frequently sampled tree is equal to the true tree: 0.79"
[1] "true tree sampled: 786 times"
```
## Whales
- Fixed tree (mean tree): ESS 10%
- Not fixed tree: ESS 1%

# Comparison of distance matrix between MrBayes trees and Bistro

## Simulates Whales
First, we need to run `mb2badger` on MrBayes output:
```
mb sim-whales-nsites1500-sorted.nex
../../Code/badger/mb2badger sim-whales-nsites1500-sorted.nex.run1.t
```
Now, we will run the `distances` executable:
```
cd Bistro/Examples/Simulations

```
- `bistro5-whales.treeBL`

## Passeriformes
- ESS fixed tree: `2.2%`
- ESS not fixed tree: `0.69%`
In this case, the gamma sampler was used 34/38 times (0.89), the truncated normal is used 3/38 times (0.08) and the exponential is used 1/38 (0.03)

## Datasets
(from `AbererStamatakisRonquist2016`)
- Run 5 datasets: `024,027,036,041,043`
  - `027` had an error with `nan` because of U instead of T, fixed!
  - `036,041` had `nan` in all branch lengths because of `?` in sequence, fixed!
  - `024,043` run without any problem, so we ran with fixed tree as well

- Results `024`:
  - ESS 0.1% (time ~9hrs)
  - ESS fixed mean tree 0.41% (time ~8hrs)
  - ESS (generalized dirichlet) 0.1%
  - ESS fixed tree (generalized dirichlet) 4.6%

- Results `027`:
  - ESS 0.13%
  - ESS fixed mean tree 8% (time ~1.5hrs)
  - ESS (generalized dirichlet) 0.1%
  - ESS fixed tree (generalized dirichlet) 2.4%

- Results `036`:
  - ESS 0.1%
  - ESS fixed mean tree 16% (time ~2.5hrs)
  - ESS (generalized dirichlet) 0.1%
  - ESS fixed mean tree (gen dir) 17%

- Results `041`:
  - ESS 0.12%
  - ESS fixed mean tree 7% (time ~2hrs)

- Results `043`:
  - ESS 0.19% (time ~ 3hrs)
  - ESS fixed mean tree 13% (time ~3hrs)
  - ESS (generalized dirichlet) 0.1%
  - ESS fixed mean (gen dir) tree 4.8%

- It is strange that the performance of bistro (on fixed tree) is so different between `024` and `043`. We checked the `.sampler` file, and in both datasets we are only using the gamma distribution to sample branch lengths (`100%`). **The reason is rates** The rates sampling is very bad in `024` (not in any of the others). We will implement the generalized dirichlet for this problem.
- Generalized Dirichlet does gives a better performance in `024`, but a little worse in others

#### After burnin = 200 & generalized dirichlet
Actually, I thought that we used 200 burnin, but I did not recompile in Birge machine, so it was using the old burnin=100
```
(master) $ grep "ESS = " *bistro*-3-*.log
birge-bistro-3-024.log:ESS = 1.95, or 0.19 percent.
birge-bistro-3-027.log:ESS = 1.00, or 0.10 percent.
birge-bistro-3-036.log:ESS = 1.00, or 0.10 percent.
birge-bistro-3-043.log:ESS = 1.48, or 0.15 percent.
birge-bistroFixT-3-024.log:ESS = 70.18, or 7.02 percent.
birge-bistroFixT-3-027.log:ESS = 85.64, or 8.56 percent.
birge-bistroFixT-3-036.log:ESS = 245.24, or 24.52 percent.
birge-bistroFixT-3-043.log:ESS = 181.98, or 18.20 percent.
bistro-3-024.log:ESS = 1.13, or 0.11 percent.
bistro-3-027.log:ESS = 1.03, or 0.10 percent.
bistro-3-036.log:ESS = 1.00, or 0.10 percent.
bistro-3-043.log:ESS = 1.85, or 0.19 percent.
bistroFixT-3-024.log:ESS = 437.87, or 43.79 percent.
bistroFixT-3-027.log:ESS = 44.28, or 4.43 percent.
bistroFixT-3-036.log:ESS = 84.44, or 8.44 percent.
bistroFixT-3-043.log:ESS = 191.86, or 19.19 percent.

(master) $ grep "ESS = " *bistro*-whales*.log
birge-bistro-whales.log:ESS = 1.55, or 0.15 percent.
birge-bistroFixT-whales.log:ESS = 62.17, or 6.22 percent.
bistro-whales.log:ESS = 2.91, or 0.29 percent.
bistroFixT-whales.log:ESS = 76.49, or 7.65 percent.

(master) $ grep "ESS = " *bistro*-cats*.log
birge-bistro-cats.log:ESS = 3.72, or 0.37 percent.
birge-bistroFixT-cats.log:ESS = 265.86, or 26.59 percent.
bistro-cats.log:ESS = 4.51, or 0.45 percent.
bistroFixT-cats.log:ESS = 46.19, or 4.62 percent.
```

Again, with burnin=200:
```
[claudia@darwin00] (10)$ grep ESS bistroFixT-4*.log
bistroFixT-4-024.log:ESS = 197.70, or 19.77 percent.
bistroFixT-4-027.log:ESS = 38.92, or 3.89 percent.
bistroFixT-4-036.log:ESS = 51.75, or 5.18 percent.
bistroFixT-4-043.log:ESS = 140.70, or 14.07 percent.
bistroFixT-4-050.log:ESS = 29.24, or 2.92 percent.
bistroFixT-4-059.log:ESS = 65.45, or 6.54 percent.
bistroFixT-4-064.log:ESS = 17.66, or 1.77 percent.
```
See the plots `convergence*.pdf` to check the convergence of loglik, pi and rates.

## Compare bistro and mb mean tree
### 024
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined024.pdf`. Distance between bistro and mb mean trees = 0.121764. On average, combined trees are a 0.0922 distance. But bistro trees are on average at 0.0429623. In fact, 97.5% of trees have distance less of 0.0544. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.061943,((((2:0.0601881,3:0.0592558):0.0295168,((4:0.0604012,5:0.0616232):0.0268623,((6:0.0609094,7:0.0564196):0.0257563,(8:0.055568,9:0.062067):0.0256383):0.0107998):0.0072168):0.00954059,(((((10:0.0584351,11:0.0481539):0.0274329,(12:0.0507255,13:0.053615):0.0202671):0.0104546,(14:0.0532989,15:0.059372):0.0261837):0.00766253,(16:0.0549529,17:0.0606936):0.0281575):0.00898972,((18:0.058076,19:0.0541883):0.0220081,(20:0.0523336,21:0.0627729):0.0288435):0.013281):0.012184):0.0145548,(22:0.0680999,23:0.053197):0.0253921):0.0232983,24:0.0542719);")
> tre.bistro = read.tree(text="(1:0.0402852,((((2:0.0396609,3:0.0392467):0.00437999,(4:0.0398325,5:0.0425959):0.00390973):0.00138546,(((6:0.0424069,7:0.0417964):0.00505505,(8:0.0411127,(9:0.0377531,10:0.0376897):0.00829165):0.00240542):0.0010178,((((11:0.0361696,12:0.037094):0.00599679,(13:0.0391003,14:0.0395072):0.00664326):0.00192441,(15:0.0437715,(16:0.0390453,17:0.0411292):0.00719562):0.00172671):0.000532715,((18:0.042361,19:0.0387008):0.0046463,(20:0.0412141,21:0.0417332):0.0058016):0.00183854):0.000604983):0.00033767):0.00159217,(22:0.0415971,23:0.0400248):0.00438032):0.00636168,24:0.0397541);")
pdf("combined024-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```

```R
> mat = as.matrix(read.table("combined024.distances"))
> mean(mat)
[1] 0.09229893
> sd(mat)
[1] 0.0506188
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.0421724 0.1243885 0.1348650 0.1565490 
> mat = as.matrix(read.table("bistro-6-024.bootstrapBL.distances"))
> mean(mat)
[1] 0.0429623
> sd(mat)
[1] 0.008153844
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
       25%        50%        75%      97.5% 
0.03897043 0.04410215 0.04853067 0.05442400
```


### 027
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined027.pdf`. Distance between bistro and mb mean trees = 0.0420587. On average, combined trees are a 0.0384 distance. But bistro trees are on average at 0.02797. In fact, 97.5% of trees have distance less of 0.0352. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.0101192,((2:0.0091657,(3:0.0115832,(4:0.00968702,(5:0.00651467,(6:0.0114666,(7:0.00812869,8:0.0235345):0.0103157):0.00762085):0.00705958):0.00458583):0.00794148):0.00422824,((((((((((((((((9:0.00918607,10:0.00851567):0.0110119,11:0.00521696):0.00640595,12:0.00973822):0.00851911,13:0.0104624):0.00839615,14:0.00242226):0.00482617,15:0.00835657):0.00905479,16:0.00967589):0.00665609,17:0.014226):0.0061813,18:0.00207056):0.00250304,19:0.0167871):0.00889612,20:0.0202873):0.00682819,21:0.0123987):0.00740288,22:0.00869153):0.00527526,23:0.0132916):0.00486893,24:0.0104183):0.0114139,(25:0.0182022,26:0.0145785):0.00584751):0.0157665):0.0069554,27:0.0128853);")
> tre.bistro = read.tree(text="(1:0.0109314,((2:0.010416,(3:0.0112083,((4:0.0106111,(5:0.00901583,6:0.0140031):0.000961532):0.00103343,7:0.0151813):0.00230838):0.00222732):0.000747447,(((((((8:0.0162405,9:0.00619896):0.000895582,10:0.00739793):0.00572367,(11:0.00915009,12:0.00821265):0.00146334):0.00286269,(13:0.0126711,(14:0.00397731,15:0.0086628):0.00223734):0.00181886):0.00190224,(((16:0.0107367,17:0.0143156):0.00172852,(18:0.00218419,19:0.0157774):0.00115909):0.00157069,(20:0.019461,21:0.012185):0.00295481):0.00076808):0.00102178,((22:0.00819718,24:0.0110085):0.00140469,23:0.0119386):0.00112501):0.0016424,(25:0.0167008,26:0.0123764):0.00164598):0.00405585):0.00203,27:0.00783759);")
pdf("combined027-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```

```R
> mat = as.matrix(read.table("combined027.distances"))
> mean(mat)
[1] 0.03839355
> sd(mat)
[1] 0.0125421
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.0273876 0.0379279 0.0501081 0.0554678
> mat = as.matrix(read.table("bistro-6-027.bootstrapBL.distances"))
> mean(mat)
[1] 0.02797086
> sd(mat)
[1] 0.004805236
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
       25%        50%        75%      97.5% 
0.02538110 0.02851375 0.03108015 0.03523185
```


### 027
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined036.pdf`. Distance between bistro and mb mean trees = 0.198598. On average, combined trees are a  distance of 0.1583. But bistro trees are on average at 0.1274. In fact, 97.5% of trees have distance less of 0.1485. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.0875644,(2:0.0834653,(3:0.0720429,(4:0.0734073,(5:0.0742511,(6:0.0839117,(7:0.0932401,(8:0.0813149,(9:0.0894596,(10:0.090185,(11:0.0876109,(12:0.0847255,(13:0.0826326,(14:0.0707776,(15:0.0781684,16:0.114987):0.0370283):0.0304587):0.0368275):0.0389148):0.029271):0.0295913):0.0208889):0.0326336):0.0368634):0.0249377):0.024963):0.0358493):0.0285766):0.0337882,(((((((((((((((((((17:0.105606,18:0.0738661):0.0363304,19:0.0765767):0.0309133,20:0.0738644):0.0384648,21:0.0821027):0.0387762,22:0.0853945):0.0308935,23:0.082245):0.0398329,24:0.0938782):0.0364868,25:0.0912975):0.0283365,26:0.09775):0.0441461,27:0.117837):0.0310755,28:0.0914752):0.0389788,29:0.0776094):0.0380291,30:0.0886566):0.0343545,31:0.0803147):0.0358186,32:0.0785755):0.0338153,33:0.084561):0.0355285,34:0.0825571):0.0093369,35:0.0721286):0.0185827,36:0.0708893):0.0275227);")
> tre.bistro = read.tree(text="(1:0.0912461,(((2:0.0834547,3:0.0776173):0.00763851,(((4:0.0757175,5:0.0748022):0.00760912,6:0.0857814):0.00230024,((7:0.0940647,(8:0.078938,9:0.0860795):0.010847):0.00213721,((10:0.0913831,11:0.0882434):0.00634471,(12:0.093142,((13:0.0863826,(14:0.0828466,15:0.0869275):0.00751384):0.00483575,(16:0.0988856,17:0.088096):0.00915454):0.00540457):0.00448426):0.00229933):0.00256046):0.0018536):0.00281701,(((((18:0.086234,19:0.0839358):0.00927216,20:0.0826623):0.00683275,(21:0.0896827,(22:0.0907858,23:0.089788):0.0149828):0.00828588):0.00890293,((((24:0.10083,(25:0.0925138,26:0.106494):0.00577601):0.00820069,27:0.120528):0.00528969,28:0.0970836):0.00311987,(29:0.0931229,30:0.0932744):0.00672053):0.00190726):0.00122616,((31:0.085475,32:0.0861242):0.00980249,33:0.0879633):0.00256983):0.00185064):0.00187825,((34:0.078133,35:0.070571):0.00786924,36:0.0708674):0.00367123);")
pdf("combined036-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```

```R
> mat = as.matrix(read.table("combined036.distances"))
> mean(mat)
[1] 0.1583718
> sd(mat)
[1] 0.06261172
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.0936739 0.1600700 0.2207330 0.2376280
> mat = as.matrix(read.table("bistro-6-036.bootstrapBL.distances"))
> mean(mat)
[1] 0.1273807
> sd(mat)
[1] 0.01788783
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.1223650 0.1302825 0.1367742 0.1485437
```

### 043
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined043.pdf`. Distance between bistro and mb mean trees = 0.101826. On average, combined trees are a 0.119 distance. But bistro trees are on average at 0.09. In fact, 97.5% of trees have distance less of 0.1039. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.0550771,((2:0.0468695,3:0.0658761):0.0195383,(((((4:0.0671356,(5:0.0401607,6:0.0634503):0.016199):0.0156356,((7:0.0527217,8:0.0368104):0.0222124,(9:0.0398978,(10:0.0262045,11:0.061646):0.0211495):0.00980919):0.00876083):0.0134813,(12:0.0533258,13:0.0477652):0.0357106):0.0102502,(((((14:0.0727,15:0.0588482):0.0231021,16:0.0657255):0.01606,(17:0.0568474,(18:0.0571245,19:0.049513):0.0229545):0.0251803):0.0121913,(((((20:0.0600766,21:0.0283899):0.0254293,22:0.0353944):0.0321892,23:0.0474676):0.0219248,24:0.0648215):0.0180102,(25:0.0416776,26:0.0507775):0.0279232):0.0175052):0.0110421,((((27:0.0407362,28:0.0157944):0.0270986,29:0.0320519):0.0208188,30:0.030959):0.0216741,((31:0.0432691,32:0.0486404):0.00879781,(33:0.0373321,34:0.0488427):0.0170062):0.0028798):0.0137634):0.00175282):0.0033213,((35:0.0569794,(36:0.0432571,37:0.0587669):0.0157889):0.0177007,(38:0.0613824,(39:0.0532691,((40:0.0593077,41:0.0743165):0.0100966,42:0.069426):0.0155246):0.00966818):0.00824363):0.0068727):0.0114034):0.0170826,43:0.0744319);")
> tre.bistro = read.tree(text="(1:0.0544328,((((((((3:0.0585844,4:0.0654338):0.00350466,((5:0.0447746,6:0.0591256):0.00533132,(7:0.0468165,8:0.0412618):0.0105817):0.00284956):0.00184366,(((9:0.0397721,10:0.0326188):0.00632344,11:0.0548882):0.00452885,(12:0.0527359,13:0.0471185):0.00812675):0.00164977):0.0023041,((14:0.0651105,(15:0.061296,16:0.0657605):0.00865914):0.00622553,(17:0.0550363,(18:0.0644016,19:0.0420601):0.0115884):0.0110939):0.00158135):0.000512852,(((((20:0.0523793,(21:0.0387686,22:0.0361629):0.0140272):0.0116357,(23:0.051239,24:0.0609308):0.00524585):0.0036179,(25:0.048484,26:0.0448579):0.00850686):0.0014545,((27:0.0271295,28:0.0272599):0.00893159,(29:0.039264,30:0.032288):0.00493417):0.0028631):0.000501311,((31:0.0411274,32:0.043376):0.00374711,33:0.0435459):0.00343855):0.000590609):0.000924653,(34:0.0421418,(35:0.0529444,(36:0.0448221,37:0.0546131):0.00982281):0.0067966):0.00107545):0.00142511,((38:0.0560059,39:0.0496388):0.00342543,((40:0.0624241,42:0.0623289):0.00292385,41:0.0709147):0.0063798):0.00228647):0.00373884,2:0.0483406):0.00456917,43:0.0653997);")
pdf("combined043-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```

```R
> mat = as.matrix(read.table("combined043.distances"))
> mean(mat)
[1] 0.1185649
> sd(mat)
[1] 0.028335
> min(mat)
[1] 0
> max(mat)
[1] 0.174224
> quantile(mat,probs=c(0.25,0.5,0.75))
      25%       50%       75% 
0.0949245 0.1268740 0.1405990 
> mat = as.matrix(read.table("bistro-6-043.bootstrapBL.distances"))
> mean(mat)
[1] 0.08638934
> sd(mat)
[1] 0.01258516
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
       25%        50%        75%      97.5% 
0.08144508 0.08776760 0.09346110 0.10386830
```

### 050
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined050.pdf`. Distance between bistro and mb mean trees = 0.049617 . On average, combined trees are a 0.0635 distance. But bistro trees are on average at 0.055. In fact, 97.5% of trees have distance less of 0.066. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.0171362,2:0.0149846,(((((((3:0.0230808,5:0.0151482):0.00842359,4:0.0119057):0.0166576,6:0.0086044):0.0118557,7:0.00960745):0.00525933,8:0.0162969):0.0152823,((((9:0.0170827,(10:0.0130696,11:0.0187387):0.0154593):0.0116027,((14:0.0160912,(15:0.0203131,16:0.0190439):0.00185829):0.00958711,((17:0.012277,(18:0.00800205,19:0.0213748):0.0122084):0.0124438,(((20:0.0215345,21:0.0157843):0.0156197,(22:0.0055496,23:0.0147809):0.0104632):0.018539,(24:0.00962791,(25:0.00764659,(26:0.0127282,(27:0.00359531,28:0.00961829):0.00857436):0.0166506):0.00519731):0.00874766):0.00565243):0.00496276):0.00380042):0.00192734,(12:0.0190936,13:0.00830736):0.0104789):0.00577314,((((((((29:0.0176279,30:0.00844828):0.00476738,31:0.0143533):0.00378477,32:0.0134895):0.00372045,(33:0.0101756,34:0.00628684):0.00428373):0.00341289,35:0.00697611):0.0104917,(((((36:0.0073245,37:0.0119651):0.00774787,(38:0.0101041,39:0.0145481):0.00594955):0.00445311,40:0.0183888):0.00544377,41:0.0120849):0.00238889,42:0.00776903):0.00254165):0.00117585,46:0.0115113):0.00154603,(((43:0.00454662,45:0.00661256):0.00367668,(47:0.00517259,48:0.010213):0.00233283):0.00455557,44:0.00478237):0.00375908):0.00514869):0.00534965):0.00382775,(49:0.0106843,50:0.0108057):0.00237488):0.00761359);")
> tre.bistro = read.tree(text="(1:0.0154804,(((((((((9:0.0175859,(10:0.0153065,11:0.0138734):0.0115353):0.00656836,(12:0.0190919,13:0.00808746):0.00453066):0.00110217,(((14:0.0151799,15:0.0214181):0.0032526,16:0.0150021):0.00194717,(17:0.0144839,(18:0.0102104,19:0.0217753):0.00732772):0.00458331):0.000699654):0.000502174,(((20:0.0240453,21:0.0147137):0.00815856,(22:0.0118928,23:0.00498894):0.00983684):0.00603102,(24:0.00731324,25:0.012913):0.00825073):0.00140745):0.00066676,(((26:0.0153121,(27:0.00722363,28:0.00796835):0.00452735):0.00422566,29:0.0173781):0.00253986,((30:0.00868715,(31:0.0147653,32:0.00901191):0.0025418):0.00137536,((33:0.00754015,34:0.00649483):0.00311134,35:0.00662817):0.000505048):0.0014893):0.00133539):0.00153713,((36:0.00520807,37:0.0134934):0.00347125,((38:0.00962268,39:0.0141955):0.00637918,40:0.0135978):0.00087323):0.000573998):0.000340206,((((41:0.00839226,42:0.00447788):0.0012356,46:0.0105008):0.000539339,((43:0.00166416,45:0.00422451):0.000669651,(47:0.00329849,48:0.00948983):0.001294):0.00117801):0.000493497,44:0.00462196):0.000324838):0.00272726,(3:0.0185381,(((4:0.0176493,5:0.0178119):0.0036971,6:0.0147107):0.00562982,(7:0.0125443,8:0.0155806):0.00376525):0.00230696):0.00431158):0.00190354,(49:0.00902211,50:0.00882045):0.00242195):0.00196616,2:0.0131855);")
pdf("combined050-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```

```R
> mat = as.matrix(read.table("combined050.distances"))
> mean(mat)
[1] 0.06359097
> sd(mat)
[1] 0.01341239
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.0544838 0.0664737 0.0737939 0.0819594 
> mat = as.matrix(read.table("bistro-6-050.bootstrapBL.distances"))
> mean(mat)
[1] 0.05483522
> sd(mat)
[1] 0.008125361
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
       25%        50%        75%      97.5% 
0.05132143 0.05558875 0.05956515 0.06620189
```

### 059
Mb and bistro trees seem to be separate. Also see the clouds of MDS space `combined059.pdf`. Distance between bistro and mb mean trees = 0.164775. On average, combined trees are a  distance of 0.184. But bistro trees are on average at 0.123. In fact, 97.5% of trees have distance less of 0.141. So, the cloud of bistro trees are more condensed among themselves than with mb tree.
```R
> tre.mb = read.tree(text="(1:0.104164,2:0.0928128,(3:0.0978276,((4:0.103273,5:0.135481):0.0207336,((6:0.107998,(7:0.121398,8:0.114914):0.0200857):0.0111001,(9:0.136844,((10:0.126164,11:0.0996231):0.0305566,(12:0.11465,(13:0.104103,((14:0.0880185,15:0.0820851):0.0330803,((16:0.0746874,17:0.0864537):0.0280768,((18:0.0759462,19:0.0800772):0.0191154,(20:0.0874527,((21:0.0808989,22:0.0815907):0.0295736,(23:0.0863543,((24:0.0744851,25:0.0800574):0.0227849,(26:0.0632752,(27:0.0719471,((28:0.0562635,29:0.0405657):0.024649,(30:0.0423431,(31:0.0283553,(32:0.0337221,(33:0.0311773,(34:0.0270126,(35:0.0422459,(36:0.0417071,(37:0.033955,(38:0.0274905,(39:0.037444,(40:0.0389946,(41:0.0364718,(42:0.0556327,(43:0.0627628,(44:0.0728088,(45:0.055417,(46:0.0648716,(47:0.0556923,(48:0.0525329,(49:0.0643463,(50:0.054136,(51:0.0482549,(52:0.0570157,(53:0.0635859,(54:0.0783192,(55:0.0691725,(56:0.0833823,(57:0.077328,(58:0.0772964,59:0.0866857):0.0235668):0.0226924):0.0209013):0.0214334):0.0238013):0.0158851):0.0236515):0.0246402):0.022559):0.02378):0.0169672):0.0215294):0.0216627):0.0126135):0.022786):0.0204875):0.0157642):0.028032):0.0146382):0.0269038):0.0252314):0.0203488):0.0243774):0.0245352):0.0221969):0.0307435):0.0292685):0.0184482):0.0200212):0.0195785):0.0187811):0.0161245):0.0155207):0.0174609):0.019148):0.0183381):0.0168219):0.0167592):0.020876):0.0194232):0.0226362):0.0138064):0.0211846):0.0165259):0.0189386):0.0264253);")
> tre.bistro = read.tree(text="(1:0.0943557,(((((5:0.121218,((6:0.10374,7:0.11369):0.00322209,(8:0.112857,9:0.118699):0.00358394):0.00304334):0.0018946,(((10:0.117187,11:0.0974169):0.00955096,(12:0.102382,(13:0.0928,(14:0.0795557,15:0.0843411):0.015092):0.00739273):0.00379112):0.00196397,((((16:0.0775687,17:0.0792895):0.00790989,(18:0.0702149,19:0.0772338):0.00374199):0.00362082,20:0.0829847):0.00397459,((((21:0.0791697,22:0.074951):0.00955584,23:0.0772221):0.00331091,((24:0.0729219,25:0.0747915):0.00410854,(26:0.057537,27:0.0664891):0.00622507):0.0014351):0.00123199,(((28:0.049463,29:0.0470755):0.0132182,((30:0.0437322,31:0.0364324):0.00997928,32:0.0468354):0.00444606):0.00238155,((33:0.0364187,34:0.0354264):0.0123958,(((35:0.0456769,36:0.043296):0.00983322,(37:0.0435059,38:0.040032):0.00797846):0.00480054,(39:0.0398097,40:0.0415728):0.0136105):0.00486356):0.00247566):0.000758158):0.00146155):0.00143918):0.000782575):0.00090563,(3:0.0957617,4:0.103788):0.00627521):0.000675562,((((41:0.0410288,42:0.0558677):0.00932555,(43:0.061924,(44:0.0664965,45:0.057146):0.0063394):0.00292283):0.00325702,(46:0.0658467,(47:0.0618802,48:0.0521427):0.00811577):0.00411472):0.00251791,((49:0.0659219,50:0.0617816):0.00741821,(51:0.0566688,(52:0.0607217,53:0.0643466):0.00676339):0.00365453):0.00348428):0.00102294):0.00141298,(((54:0.0764412,55:0.0699061):0.00729012,(56:0.081546,57:0.0755053):0.0105512):0.00414401,(58:0.0794992,59:0.0787373):0.0105822):0.00274449):0.00625678,2:0.0965757);")
pdf("combined059-trees.pdf")
layout(m=matrix(c(1,2),nrow=1))
plot(tre.mb, no.margin=TRUE)
plot(tre.bistro, no.margin=TRUE)
dev.off()
```
The Mb mean tree is very weird, let's compare to the MAP:
```R
tre.map = read.tree(text="(2,((((((((((((((((((((((((((((((((((((((((((((((59,58),57),56),55),54),53),52),51),50),49),48),47),46),45),44),43),42),41),40),39),38),37),36),35),34),33),32),31),30),(29,28)),27),26),(25,24)),23),(22,21)),20),(19,18)),(17,16)),(15,14)),13),12),(11,10)),9),((8,7),6)),(5,4)),3),1);")
```
Same!

```R
> mat = as.matrix(read.table("combined059.distances"))
> mean(mat)
[1] 0.1840661
> sd(mat)
[1] 0.06415083
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
     25%      50%      75%    97.5% 
0.122241 0.220355 0.231730 0.284691 
> mat = as.matrix(read.table("bistro-6-059.bootstrapBL.distances"))
> mean(mat)
[1] 0.1232917
> sd(mat)
[1] 0.01570512
> quantile(mat,probs=c(0.25,0.5,0.75,0.975))
      25%       50%       75%     97.5% 
0.1191770 0.1252930 0.1308088 0.1408361
```

### Generalized Dirichlet

I think the natural thing to do is to use

X_i ~ Gamma(alpha_i, lambda_i)  (Note lambda_i may not equal lambda_j).

Let Y_i  = X_i / S  where S = sum(X_i)

There is a literature on the distribution of S in this general case. Maybe there is one for the vector of Y_i too.

I am not sure how to fit this, but if we cannot make the mean scale work, this is the direction I think we should try.

-Bret

------------------------
I have an expression for the density of x_1,...,x_k when

x_i = y_i/sum(y_j)

and y_i ~ iid Gamma(alpha_i,lambda_i).

So, we have a density and can generate random variables from it. But I am getting stuck trying to find a closed expression for the mean (and hence the variance too), so I am not quite sure how we could find good alpha and lambda values to try to match the empirical values.

As described, the parameters are not identifiable, but requiring the lambda parameters to have mean 1 (or maybe sum to one?) would make it work.

-Bret

------------------------------
I was able to use the mean and variances from the MCMC for the s from the 024 data set. I used these to calculate alpha and lambda for each component. I then generated random gammas using these variables and renormalized for proposed s values. This is a generalized Dirichlet, but not like the literature. Not sure if it is new or not, but I have not found it in the literature. The triplots for these generated s values match the ones from MCMC very well.

I will try to code this for an alternative approach and see if this helps.

-Bret


### Inf in 041

I found a problem in `041` where we have `logBL=inf`:
```shell
grep -n "inf" *
bistroFixT041---0-249.out:101:(1,(((((2,5),(9,24)),39),((((((3,10),(4,(7,41))),((6,25),(18,26))),(((11,15),23),(33,40))),((8,(31,32)),(13,22))),(((((12,20),(29,30)),(((19,28),27),21)),37),((14,36),(34,35))))),(17,38)),16); -16228.6  0 inf 150.701 21.8951 -inf 0.308224 0.183829 0.227101 0.280846 0.0957594  0.206061  0.139481  0.107697  0.311803  0.139198
bistroFixT041---250-499.out:240:(1,(((((2,5),(9,24)),39),((((((3,10),(4,(7,41))),((6,25),(18,26))),(((11,15),23),(33,40))),((8,(31,32)),(13,22))),(((((12,20),(29,30)),(((19,28),27),21)),37),((14,36),(34,35))))),(17,38)),16); -16217.8  0 inf 150.877 29.4758 -inf  0.30887  0.17697 0.235047 0.279112 0.0992608  0.223081  0.150875 0.0850139  0.309227  0.132542
bistroFixT041---500-749.out:39:(1,(((((2,5),(9,24)),39),((((((3,10),(4,(7,41))),((6,25),(18,26))),(((11,15),23),(33,40))),((8,(31,32)),(13,22))),(((((12,20),(29,30)),(((19,28),27),21)),37),((14,36),(34,35))))),(17,38)),16); -16222.5  0 inf 151.134 26.8807 -inf 0.310801 0.184869 0.248141  0.25619 0.108572 0.202526 0.144592 0.088879 0.314319 0.141112
bistroFixT041---500-749.out:44:(1,(((((2,5),(9,24)),39),((((((3,10),(4,(7,41))),((6,25),(18,26))),(((11,15),23),(33,40))),((8,(31,32)),(13,22))),(((((12,20),(29,30)),(((19,28),27),21)),37),((14,36),(34,35))))),(17,38)),16); -16218.8  0 inf 151.219 26.6197 -inf 0.305372 0.175671 0.235902 0.283056  0.118523  0.197061  0.162325 0.0876721  0.316597  0.117822
```
But if we look in the `treeBL`, the branch lengths do not seem to be weird in any sense!
```
((39:0.0371741,((((37:0.0653909,(((12:0.0431819,20:0.0689953):0.0164853,(30:0.0723324,29:0.1372519):0.0324667):0.0034364,(((19:0.0916155,28:0.0554377):0.0197495,27:0.0488757):0.0176752,21:0.0907134):0.0078632):0.0165991):0.0020457,((34:0.0900311,35:0.0718004):0.0333207,(14:0.0338025,36:0.0521658):0.0171395):0.0025286):0.0099271,((((40:0.0648464,33:0.0695947):0.0196077,((11:0.0790237,15:0.0635204):0.0279725,23:0.0656730):0.0158994):0.0008214,(((3:0.0490240,10:0.0663798):0.0388326,(4:0.1115355,(41:0.0784227,7:0.0525318):0.0149854):0.0239091):0.0038172,((26:0.0223549,18:0.0431359):0.0183500,(6:0.0966291,25:0.0469873):0.0106713):0.0167755):0.0039002):0.0021800,((22:0.0838966,13:0.0358311):0.0186856,(8:0.0762191,(32:0.0770553,31:0.0482529):0.0022179):0.0155292):0.0008749):0.0011405):0.0039442,((38:0.0530934,17:0.0359830):0.0109768,(16:0.0572615,1:0.0697440):0.0298572):0.0076620):0.0000552):0.0025057,(9:0.0707964,24:0.0432450):0.0072144,(2:0.0693594,5:0.0714157):0.0027924);
```

## Distance matrix (plot `combined.pdf`, script `mds2.R`)

- Artiodactyl: MB and Bistro cloud are intersected, mean trees are close.
- Cats and dogs: separate clouds for different topologies, MB and Bistro separate clouds, distant means
- Whales: MB and bistro different clouds, not separated by topology, distant means
- 043: MB bimodal, different cloud for bistro, distant means
- 024: Bisto does not agree with any mb tree, MB trimodal, distant means

**Conclusion** The pull towards the mean tree does not work because bistro mean tree tends to be different (and far) from MB mean tree.

## Improvements to MCMC

We identified issues with MCMC in which the specified burnin was not enough. We needed some measures inside the code to determine whether the chain has converged.
We modified the mcmc step to do 4 independent chains, and calculate the Gelman-Rubin statistics to decide if the chains have converged (using 50% of the chain as burnin).

Performance on fixed tree:
- cats dogs: ESS = 119.57, or 11.96 percent.
- whales: ESS = 36.60, or 3.66 percent.
```
[claudia@darwin03] (23)$ grep ESS bistroFixT-6*.log
bistroFixT-6-024.log:ESS = 381.25, or 38.12 percent.
bistroFixT-6-027.log:ESS = 32.99, or 3.30 percent.
bistroFixT-6-036.log:ESS = 55.41, or 5.54 percent.
bistroFixT-6-043.log:ESS = 110.77, or 11.08 percent.
bistroFixT-6-050.log:ESS = 3.91, or 0.39 percent.
bistroFixT-6-059.log:ESS = 166.70, or 16.67 percent.
```

Before:
```
[claudia@darwin00] (10)$ grep ESS bistroFixT-4*.log
bistroFixT-4-024.log:ESS = 197.70, or 19.77 percent.
bistroFixT-4-027.log:ESS = 38.92, or 3.89 percent.
bistroFixT-4-036.log:ESS = 51.75, or 5.18 percent.
bistroFixT-4-043.log:ESS = 140.70, or 14.07 percent.
bistroFixT-4-050.log:ESS = 29.24, or 2.92 percent.
bistroFixT-4-059.log:ESS = 65.45, or 6.54 percent.
bistroFixT-4-064.log:ESS = 17.66, or 1.77 percent.
```

## Mixture of bootstrap and MCMC samples
Not a great performace:
```
[claudia@darwin03] (8)$ grep "ESS" bistro-7*.log
bistro-7-024.log:ESS = 1.00, or 0.10 percent.
bistro-7-027.log:ESS = 1.00, or 0.10 percent.
bistro-7-036.log:ESS = 1.00, or 0.10 percent.
bistro-7-043.log:ESS = 1.00, or 0.10 percent.
bistro-7-050.log:ESS = 1.94, or 0.19 percent.
bistro-7-059.log:ESS = 1.00, or 0.10 percent.
```
We will study if the MB tree has a big weight in the sample, or not.
### Whales
We have two MB trees, with only difference in clades:
- `((12,10),(13,(11,(9,8))))` with PP = 0.15
- `((13,(12,11)),(10,(9,8)))` with PP = 0.11

The tree with minimum parsimony in `bistro4.pmap`:
```
(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((((16,17),((18,19),(20,21))),31),(((24,(27,28)),(25,26)),(29,30))),(22,23))))); 0.917937
```
which is very close to the second tree in MB, but not identical.
```R
library(ape)
mbtre1 = read.tree(text="(2,(((((31,((30,29),((26,25),((28,27),24)))),(23,22)),(((21,20),(19,18)),(17,16))),15),(((((12,10),(13,(11,(9,8)))),(7,6)),(14,(5,4))),3)),1);") ## PP = 0.15
mbtre2 = read.tree(text="(2,(((((31,((30,29),((26,25),((28,27),24)))),(23,22)),(((21,20),(19,18)),(17,16))),15),(((((13,(12,11)),(10,(9,8))),(7,6)),(14,(5,4))),3)),1);") ## PP = 0.11
minPars = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((((16,17),((18,19),(20,21))),31),(((24,(27,28)),(25,26)),(29,30))),(22,23)))));")
layout(matrix(c(1,2),nrow=1))
plot(mbtre2, no.margin=TRUE)
plot(minPars, no.margin=TRUE)
```
We want to know the weight that we give to the MB trees. We search for the canonical form of the MrBayes trees in `whales.nex.run1.top`:
- MB tree 1, with PP=0.15 in canonical form, does not appear in `bistro4.pmap` (from MCMC sample) nor in `bistro4.dmap` (from bootstrap sample): `(1,2,((3,(((4,5),14),((6,7),((((8,9),11),13),(10,12))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));`
- MB tree 2, with PP=0.11 in canonical form, does not appear in `bistro4.pmap` (from MCMC sample) nor in `bistro4.dmap` (from bootstrap sample): `(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));`

So, the clade `{22-31}` has high support in MB, but not in Bistro. It is not that bad in the combined sample `bistro4-dist.smap`: `0.10155057 {22-31}`, but in MCMC it has zero: `0.00000000 {22-31}`.

Still, the most probable tree in Bistro is very close to the best tree in MrBayes:
```R
bistrotre = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),(((24,((25,26),(27,28))),(29,30)),31))))));")
plot(mbtre1,no.margin=TRUE)
plot(bistrotre,no.margin=TRUE)
```
When we look at the trees with bigger weight:
```R
source("../../Scripts/readBistro.R")
dat = readBistro("bistro4")
dat[which(dat$w>0.1),]
tree
116 (1,2,((3,(((4,5),14),((6,7),((((8,9),11),13),(10,12))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));
533 (1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),(((24,((25,26),(27,28))),(29,30)),31))))));
959 (1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));
logl   logTop   logBL logPrior    logQ    logWt      pi1      pi2
116 -12700.1 -8.19213 230.458  114.071 32.0449 -12840.3 0.236586 0.224867
533 -12715.8 -9.33905 218.178  113.526 28.0682 -12839.2 0.239776 0.215767
959 -12712.5 -4.93804 215.264  113.486 30.1094 -12839.4 0.233628 0.219768
pi3      pi4       s1       s2        s3         s4       s5
116 0.174140 0.364407 0.156044 0.215387 0.0539366 0.01190150 0.554446
533 0.171681 0.372776 0.150807 0.192966 0.0647123 0.00958182 0.571554
959 0.168277 0.378328 0.160774 0.181724 0.0552653 0.01111810 0.584949
s6         w
116 0.00828461 0.1147488
533 0.01037850 0.3447246
959 0.00616993 0.2822366
```
The second tree in bistro with weight = 0.28 is almost identical to the best tree in MrBayes
```R
bistro2 = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));")
```
For this dataset, we have ESS=4% with fixed tree and ESS=0.4% for random tree.
What is special about the bistro trees? Nothing special about the logl, nor logTop, nor logPrior, nor logQ, nor logBL
```R
ggplot(dat,aes(x=logl,y=logTop,color=w))+geom_point()
ggplot(dat,aes(x=logPrior,y=logQ,color=w))+geom_point()
```

### 024

```R
library(ape)
mbtre1 = read.tree(text="(24,((23,22),((((21,20),(19,18)),((17,16),((15,14),((13,12),(11,10))))),((((9,8),(7,6)),(5,4)),(2,3)))),1);") ## PP = 0.34
mbtre2 = read.tree(text="(2,((24,23),((((((((22,21),(20,19)),(18,17)),(16,15)),((14,13),(12,11))),(10,9)),((8,7),(6,5))),(4,3))),1);") ## PP = 0.30
minPars = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((((16,17),((18,19),(20,21))),31),(((24,(27,28)),(25,26)),(29,30))),(22,23)))));")
layout(matrix(c(1,2),nrow=1))
plot(mbtre1, no.margin=TRUE)
plot(mbtre2, no.margin=TRUE)
```
These trees are very weird: it seems like if 2 was inserted int a different place, and all the other taxa was simply shifted.

**aqui voy**
- need to compare to min parsimony tree in pmap, and the remaining:
We want to know the weight that we give to the MB trees. We search for the canonical form of the MrBayes trees in `whales.nex.run1.top`:
- MB tree 1, with PP=0.15 in canonical form, does not appear in `bistro4.pmap` (from MCMC sample) nor in `bistro4.dmap` (from bootstrap sample): `(1,2,((3,(((4,5),14),((6,7),((((8,9),11),13),(10,12))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));`
- MB tree 2, with PP=0.11 in canonical form, does not appear in `bistro4.pmap` (from MCMC sample) nor in `bistro4.dmap` (from bootstrap sample): `(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));`

So, the clade `{22-31}` has high support in MB, but not in Bistro. It is not that bad in the combined sample `bistro4-dist.smap`: `0.10155057 {22-31}`, but in MCMC it has zero: `0.00000000 {22-31}`.

Still, the most probable tree in Bistro is very close to the best tree in MrBayes:
```R
bistrotre = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),(((24,((25,26),(27,28))),(29,30)),31))))));")
plot(mbtre1,no.margin=TRUE)
plot(bistrotre,no.margin=TRUE)
```
When we look at the trees with bigger weight:
```R
source("../../Scripts/readBistro.R")
dat = readBistro("bistro4")
dat[which(dat$w>0.1),]
tree
116 (1,2,((3,(((4,5),14),((6,7),((((8,9),11),13),(10,12))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));
533 (1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),(((24,((25,26),(27,28))),(29,30)),31))))));
959 (1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));
logl   logTop   logBL logPrior    logQ    logWt      pi1      pi2
116 -12700.1 -8.19213 230.458  114.071 32.0449 -12840.3 0.236586 0.224867
533 -12715.8 -9.33905 218.178  113.526 28.0682 -12839.2 0.239776 0.215767
959 -12712.5 -4.93804 215.264  113.486 30.1094 -12839.4 0.233628 0.219768
pi3      pi4       s1       s2        s3         s4       s5
116 0.174140 0.364407 0.156044 0.215387 0.0539366 0.01190150 0.554446
533 0.171681 0.372776 0.150807 0.192966 0.0647123 0.00958182 0.571554
959 0.168277 0.378328 0.160774 0.181724 0.0552653 0.01111810 0.584949
s6         w
116 0.00828461 0.1147488
533 0.01037850 0.3447246
959 0.00616993 0.2822366
```
The second tree in bistro with weight = 0.28 is almost identical to the best tree in MrBayes
```R
bistro2 = read.tree(text="(1,2,((3,(((4,5),14),((6,7),(((8,9),10),((11,12),13))))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));")
```
For this dataset, we have ESS=4% with fixed tree and ESS=0.4% for random tree.
What is special about the bistro trees? Nothing special about the logl, nor logTop, nor logPrior, nor logQ, nor logBL
```R
ggplot(dat,aes(x=logl,y=logTop,color=w))+geom_point()
ggplot(dat,aes(x=logPrior,y=logQ,color=w))+geom_point()
```
We can also see the MDS plot, and we can see that the bistro cloud intersects slightly with the MB cloud, but not enough
