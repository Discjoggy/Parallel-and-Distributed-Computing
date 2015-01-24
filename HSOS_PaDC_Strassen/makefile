CC = icc
#CC = g++
CFLAGS = -O3 -fmessage-length=0 -msse4.2 -march=native -ffast-math -fforce-addr
LDFLAGS = -ltbb -ltbbmalloc

all: HSOS_PaDC_Strassen

Definitions.o: Definitions.cpp Definitions.h
	${CC} ${CFLAGS} -c Definitions.cpp

Strassen.o: Strassen.cpp Strassen.h Matrix.h
	${CC} ${CFLAGS} -c Strassen.cpp

HSOS_PaDC_Strassen: Definitions.o Strassen.o Main.o
	${CC} ${CFLAGS} Definitions.o Strassen.o Main.o ${LDFLAGS} -o HSOS_PaDC_Strassen

clean:
	rm -rf *.o HSOS_PaDC_Strassen
