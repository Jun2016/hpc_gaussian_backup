#CC = gcc
#CCFLAGS = -fopenmp -lm

MPICC = mpicc
MPICCFLAGS = -lm

LIB_DIR = library
INCLUDE = -I$(LIB_DIR)

SRC_DIRS = step_2 step_3 step_6
SRC = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.c))
EXEC = $(SRC:.c=.x)

LIB_C = $(wildcard $(LIB_DIR)/*.c)
LIB_H = $(wildcard $(LIB_DIR)/*.h)
LIB_O = $(LIB_C:.c=.o)

#OPENMP_DIR = step_5

.PHONY: all clean

.SECONDARY:

all: $(EXEC)

$(LIB_DIR)/%.o: $(LIB_DIR)/%.c $(LIB_DIR)/%.h
	$(MPICC) -c -o $@ $< $(MPICCFLAGS)

#$(OPENMP_DIR)/%.o: $(OPENMP_DIR)/%.c ../$(LIB_H)
#	$(CC) -c -o $@ $< $(CCFLAGS)

#$(OPENMP_DIR)/%.x: $(OPENMP_DIR)/%.o ../$(LIB_H)
#	$(CC) $(CCFLAGS) -o $@ $^ -lm

%.x: $(LIB_O) %.o
	$(MPICC) $(MPICCFLAGS) -o $@ $^ -lm

%.o: %.c $(LIB_H)
	$(MPICC) -c $(INCLUDE) $(MPICCFLAGS) -o $@ $<



clean:
	find . -type f -name *.o -print -delete
	find . -type f -name *.x -print -delete
#	find . -type f -name *.mm -print -delete

