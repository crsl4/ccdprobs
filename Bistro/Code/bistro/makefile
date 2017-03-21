#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#
##OBJS = parameter.o sequence.o ccdprobs.o tree.o model.o random.o CladeTable.o SplitTable.o Tree-summarize.o Trees.o Prime.o Clade.o
OBJS = parameter.o sequence.o ccdprobs.o tree.o model.o random.o cladegraph.o
TARGET = all

#CC = g++-5
#CFLAGS = -O3 -I. -std=c++11 -Wno-deprecated

#CC = g++
#CFLAGS = -O3 -I. -std=c++11

# current:
CC = g++
CFLAGS = -O3 -I. -std=c++11 -Wno-deprecated -Wignored-attributes -Wdeprecated-declarations
LDLIBS = -lpthread

# darwin:
#CC = /s/gcc-5.3.0/bin/g++
#CFLAGS = -O3 -I. -std=c++11 -Wno-deprecated -Wdeprecated-declarations
#LDLIBS = -lpthread

# to debug with gdb, takes hours!
# CC = g++
# CFLAGS = -I. -std=c++11 -Wno-deprecated -Wignored-attributes -Wdeprecated-declarations -g -v -da -Q
# LDLIBS = -lpthread

.SUFFIXES: .C .o

all: bistro bmcmc #bccd

##bistro: main-threads.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h CladeTable.h SplitTable.h Tree-summarize.h Trees.h Prime.h Clade.h $(OBJS)
bistro: main-threads.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h cladegraph.h $(OBJS)
	$(CC) $(CFLAGS) -o bistro main-threads.C $(OBJS) $(LDLIBS)

 bmcmc: main-mcmc.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h $(OBJS)
	$(CC) $(CFLAGS) -o bmcmc main-mcmc.C $(OBJS) $(LDLIBS)

# bccd: main-ccd.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h $(OBJS)
# 	$(CC) $(CFLAGS) -o bccd main-ccd.C $(OBJS) $(LDLIBS)

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(OBJS) bistro bmcmc #bccd
