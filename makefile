tspmain:tspga.o point.o tspmain.o
	g++ tspga.o point.o tspmain.o -o tspmain

tspmain.o:tsp-main.cc
	g++ tsp-main.cc -c -Wall -g -o tspmain.o
point.o:Point.cc
	g++ Point.cc -c -Wall -g -o point.o
tspga.o:tsp-ga.cc
	g++ tsp-ga.cc -c -Wall -g -o tspga.o