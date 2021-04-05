CC = gcc
CC2 = icc
CFLAGS2 = -xW -ipo -O3 -unsafe-math-optimization -ftz
CFLAGS3 = -O3 -ffast-math -march=opteron -ftz
CFLAGS = -pg -g
LIBS = -lm
OBJS =
all: icc
main: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} multigrid.c ${LIBS}
icc: ${OBJS}
	${CC2} ${CFLAGS2} ${OBJS} multigrid.c ${LIBS}
clean:
	rm -f *~ *.o
