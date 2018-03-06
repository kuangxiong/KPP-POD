include make.inc
main:
	${MAKE} main.o
	(cd src; make all)
	$(CC) -o main ./obj/*.o $(LIBS)

.PHONY:clean
clean:
	rm -rf *.out *sig fort.* obj/*.o main
