#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#
OBJS = parameter.o sequence.o ccdprobs.o tree.o model.o random.o
TARGET = all

#CC = g++-5
#CFLAGS = -O3 -I. -std=c++11 -Wno-deprecated

#CC = g++
#CFLAGS = -O3 -I. -std=c++11

CC = g++
CFLAGS = -O3 -I. -std=c++11 -Wno-deprecated -Wignored-attributes -Wdeprecated-declarations #-g -O0

.SUFFIXES: .C .o

all: bistro

bistro: main-threads.C parameter.h sequence.h alias.h ccdprobs.h model.h tree.h random.h $(OBJS)
	$(CC) $(CFLAGS) -o bistro main-threads.C $(OBJS)

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(OBJS) bistro
