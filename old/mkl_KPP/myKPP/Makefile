include make.inc

main:
	${MAKE} main.o
	$(MAKE) -C src 
	$(CC) -o main ./obj/*.o $(LIBS) 

.PHONY:clean
clean:
	rm -rf obj/*.o main *.text
