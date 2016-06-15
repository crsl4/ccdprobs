#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#
OBJS = parameter.o sequence.o ccdprobs.o tree.o model.o
TARGET = all

#CC = g++-5
CC = g++
CFLAGS = -O3 -I. -std=c++11 #-Wno-deprecated
#CFLAGS = -I. -g

.SUFFIXES: .C .o

all: bistro

bistro: main.C parameter.h sequence.h alias.h ccdprobs.h tree.h model.h $(OBJS)
	$(CC) $(CFLAGS) -o bistro main.C $(OBJS)

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(OBJS) bistro
