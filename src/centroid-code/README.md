Downloaded from here: http://comet.lehman.cuny.edu/owen/code.html
(Megan Owen website)

If you download the "SturmMean" code, then you can run it really easy. The code needs trees with branch lengths, so I modified bistro to print the list of bootstrap trees with branch lengths: see bistro7.bootstrapBL (new bistro pushed to github, so you can do git pull):

So, in Bistro/Examples/Artiodactyl, you can run
```shell
java -jar ../../../../SturmMean/SturmMean.jar -o bistro7.sturmmean bistro7.bootstrapBL
```
(assuming that you put the SturmMean folder outside of ccdprobs), and you get the mean tree in the output file `bistro7.stummean` (which actually matches the MAP tree from MrBayes, I think).

Megan told me that we want the number of iterations to be close to ~10,000 (which in this case it is not), so we can add the option:
`-c 10` (which I don't understand well).

Anyway, from few examples that I've done, it seems that the mean tree is close (or equal) to the MLE tree. Megan suspects this is true, but no one has proven this.

Megan does not have a clean version of the code to calculate distances, she just sent me the jar file (in the folder now).
