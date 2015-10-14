# ======================
# Opcions d'optimització
# ======================
COMPILER=gcc
OPT=-g -Wall -std=c99
#OPT=-O3 -Wall -pedantic -std=c99 

all : wavefunction trajectory

# =========
# Utilitats
# =========
wavefunction : wavefunction.o numerical_lib.o quantum_lib.o
	$(COMPILER) -o wavefunction $(OPT) wavefunction.o numerical_lib.o quantum_lib.o -lm
trajectory : trajectory.o numerical_lib.o quantum_lib.o
	$(COMPILER) -o trajectory $(OPT) trajectory.o numerical_lib.o quantum_lib.o -lm

# ===========
# Bibliotecas
# ===========
wavefunction.o : wavefunction.c
	$(COMPILER) -c $(OPT) wavefunction.c -lm
trajectory.o : trajectory.c
	$(COMPILER) -c $(OPT) trajectory.c -lm
numerical_lib.o : numerical_lib.c
	$(COMPILER) -c $(OPT) numerical_lib.c -lm
quantum_lib.o : quantum_lib.c
	$(COMPILER) -c $(OPT) quantum_lib.c -lm

# ======
# Neteja
# ======
clean :
	rm -f *.o
realclean : clean
	rm -f wavefunction trajectory


