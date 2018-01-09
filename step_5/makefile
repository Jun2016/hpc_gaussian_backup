CC = gcc
CCFLAGS = -fopenmp -lm

all: 
	$(CC) $(CCFLAGS) -o step_5.x step_5.c ../library/mmio.c ../library/helper.c

clean:
	find . -name *.o -print -delete
	find . -name *.x -print -delete


