#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#

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

all: bistro distances

bistro: main-threads.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h cladegraph.h $(OBJS)
	$(CC) $(CFLAGS) -o bistro main-threads.C $(OBJS) $(LDLIBS)

#bmcmc: main-mcmc.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h $(OBJS)
#	$(CC) $(CFLAGS) -o bmcmc main-mcmc.C $(OBJS) $(LDLIBS)

# bccd: main-ccd.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h $(OBJS)
# 	$(CC) $(CFLAGS) -o bccd main-ccd.C $(OBJS) $(LDLIBS)

distances: main-distances.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h cladegraph.h $(OBJS)
	$(CC) $(CFLAGS) -o distances main-distances.C $(OBJS) $(LDLIBS)

clean:
	-rm -f $(OBJS) bistro distances

ccdprobs.o: ccdprobs.C alias.h ccdprobs.h
	$(CC) $(CFLAGS) -c ccdprobs.C

cladegraph.o: cladegraph.C cladegraph.h tree.h
	$(CC) $(CFLAGS) -c cladegraph.C

model.o: model.C sequence.h parameter.h model.h tree.h
	$(CC) $(CFLAGS) -c model.C

parameter.o: parameter.C parameter.h
	$(CC) $(CFLAGS) -c parameter.C

random.o: random.C random.h
	$(CC) $(CFLAGS) -c random.C

sequence.o: sequence.C model.h sequence.h
	$(CC) $(CFLAGS) -c sequence.C

tree.o: tree.C model.h random.h tree.h
	$(CC) $(CFLAGS) -c tree.C
