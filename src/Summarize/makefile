CC = g++

CCFLAGS = -O3
LFLAGS = 

OBJS = CladeTable.o SplitTable.o Tree.o Trees.o Prime.o Clade.o summarize.o 

all    :  summarize

summarize : ${OBJS}
	${CC} -o $@ ${CCFLAGS} ${OBJS} ${LFLAGS}

clean	: 
	rm *.o summarize
