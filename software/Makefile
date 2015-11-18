# ======================
# Opcions d'optimització
# ======================
COMPILER=gcc
#OPT=-g -Wall -std=c99
OPT=-O3 -Wall -std=c99 
TAIL=-lm

all : wavefunction trajectory area

# =========
# Utilitats
# =========
wavefunction : wavefunction.o numerical_lib.o quantum_lib.o functions.o
	$(COMPILER) -o wavefunction $(OPT) wavefunction.o numerical_lib.o quantum_lib.o functions.o $(TAIL)
trajectory : trajectory.o numerical_lib.o quantum_lib.o functions.o
	$(COMPILER) -o trajectory $(OPT) trajectory.o numerical_lib.o quantum_lib.o functions.o $(TAIL)
area : area.o numerical_lib.o quantum_lib.o functions.o
	$(COMPILER) -o area $(OPT) area.o numerical_lib.o quantum_lib.o functions.o $(TAIL)

# ===========
# Bibliotecas
# ===========
wavefunction.o : wavefunction.c
	$(COMPILER) -c $(OPT) wavefunction.c $(TAIL) 
trajectory.o : trajectory.c
	$(COMPILER) -c $(OPT) trajectory.c $(TAIL)
numerical_lib.o : numerical_lib.c
	$(COMPILER) -c $(OPT) numerical_lib.c $(TAIL)
quantum_lib.o : quantum_lib.c
	$(COMPILER) -c $(OPT) quantum_lib.c $(TAIL)
functions.o : functions.c
	$(COMPILER) -c $(OPT) functions.c $(TAIL)

# ======
# Neteja
# ======
clean :
	rm -f *.o
realclean : clean
	rm -f wavefunction trajectory area


