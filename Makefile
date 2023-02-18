# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 


CC = gcc -g
EXEC = cgSolver
CFLAG =-std=c99 -lm -O3 -mavx -march=native -DLIKWID_PERFMON -I/home/soft/likwid/include  
LFLAGS = -lm -L${LIKWID_HOME}/lib -llikwid

MODULOS = sislin \
	utils \
	resolvedorGradConjug 
 
OBJETOS = main.o $(addsuffix .o,$(MODULOS))

.PHONY: all clean

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

all: ${EXEC}
${EXEC}: ${EXEC}.o
${EXEC}: ${OBJETOS}
	$(CC) ${CFLAGS} -O $@ $^ $(LFLAGS)


clean:
	$(RM) $(OBJETOS) $(EXEC)