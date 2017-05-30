OBJS=algos.o allocate.o control1.o control2.o dofunct.o eignv.o flops.o function.o gettok.o itemcopy.o main.o matrix.o mem.o print.o qr.o stack.o twobytwo.o
matrixcalc: ${OBJS}
	gcc -o $@ -lm ${OBJS}
clean:
	rm -f ${OBJS} matrixcalc
