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
