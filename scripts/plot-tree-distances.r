foo = read.table("foo", header=true)
require(dplyr)
foo$topology = with(foo, reorder(topology,distance))
foo %>% group_by(topology) %>% summarize(mean=mean(distance),n=n())
