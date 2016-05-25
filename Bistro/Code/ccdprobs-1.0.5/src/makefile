# Bret Larget
# July 17, 2012
#
# Copyright 2012 (C)
#
# To remove files, type "make clean"
#

TARGET = all

CC = g++
CFLAGS = -O3 -I.

.SUFFIXES: .C .o 

all: ccdprobs mbsum distance

ccdprobs: ccdprobs.C
	$(CC) $(CFLAGS) -o ccdprobs ccdprobs.C

mbsum: mbsum.C mbsumtree.o
	$(CC) $(CFLAGS) -o mbsum mbsum.C mbsumtree.o

distance: distance.C distance.h tree.h
	$(CC) $(CFLAGS) -o distance distance.C

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f *.o ccdprobs mbsum distance
