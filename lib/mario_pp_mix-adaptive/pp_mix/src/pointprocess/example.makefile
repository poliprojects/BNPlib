example : example.o point2pattern.o pointprocess.o strauss.o sampler.o
	g++ -std=c++11 -o example example.o point2pattern.o pointprocess.o strauss.o sampler.o -lm -lgsl -lgslcblas -O3

point2pattern.o : point2pattern.c
	g++ -I/usr/include/eigen3 -c point2pattern.c -g

pointprocess.o : pointprocess.c
	g++ -I/usr/include/eigen3 -c pointprocess.c -g

sampler.o : sampler.c
	g++ -I/usr/include/eigen3 -c sampler.c -g

strauss.o : strauss.c
	g++ -I/usr/include/eigen3 -c strauss.c -g

example.o : example.c
	g++ -I/usr/include/eigen3 -c example.c -g

clean:
	rm *.o 
tar:
	tar -zcvf c.tar.gz *.c *.h makefile
