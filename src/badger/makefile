#CC = condor_compile g++
CC = g++

CCFLAGS = -O3 
LFLAGS = 

OBJS = util.o runsettings.o genedata.o permutation.o event.o tree.o files.o main.o gpl.o

all    : badger summarize genrc chart

badger : ${OBJS}
	${CC} -o $@ ${CCFLAGS} ${OBJS}

clean  :
	rm *.o badger summarize genrc chart

main.o : util.h runsettings.h genedata.h tree.h files.h main.C rand.h event.h
	${CC} -c ${CCFLAGS} main.C

tree.o : util.h runsettings.h genedata.h tree.h tree.C rand.h permutation.h event.h
	${CC} -c ${CCFLAGS} tree.C

genedata.o : util.h runsettings.h genedata.h genedata.C rand.h permutation.h 
	${CC} -c ${CCFLAGS} genedata.C 

permutation.o : util.h permutation.h permutation.C stack.h extent.h
	${CC} -c ${CCFLAGS} permutation.C 

files.o : runsettings.h files.h files.C
	${CC} -c ${CCFLAGS} files.C	

runsettings.o: util.h runsettings.h runsettings.C
	${CC} -c ${CCFLAGS} runsettings.C

util.o: util.h util.C
	${CC} -c ${CCFLAGS} util.C

event.o: event.h util.h rand.h permutation.h runsettings.h event.C
	${CC} -c ${CCFLAGS} event.C

gpl.o:	gpl.h gpl.C
	${CC} -c ${CCFLAGS} gpl.C

summarize: summarize.C prime.h votelist.h
	${CC} -o $@ ${CCFLAGS} summarize.C

genrc:	genrc.C
	${CC} -o $@ ${CCFLAGS} genrc.C

chart:	chart.C
	${CC} -o $@ ${CCFLAGS} chart.C
