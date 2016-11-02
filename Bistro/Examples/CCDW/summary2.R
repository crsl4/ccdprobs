## just want to identify one odd point

source("../../Scripts/readBistro.r")
bistro = readBistro("test11-3")
data = readDataSort("test11-3")
head(bistro)
bistro[bistro$w>0.6,] ## 576
data[576,]
summary(bistro)
