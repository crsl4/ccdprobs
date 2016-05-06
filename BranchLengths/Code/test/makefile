#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#
OBJS = parameter.o sequence.o random.o tree.o model.o
TARGET = all

CC = g++
CFLAGS = -O3 -I.
#CFLAGS = -I. -g

.SUFFIXES: .C .o 

all: bl

bl: main.C parameter.h sequence.h random.h tree.h model.h $(OBJS)
	$(CC) $(CFLAGS) -o bl main.C $(OBJS)

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(OBJS) bl
