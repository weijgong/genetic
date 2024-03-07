plan: test.o SAT_VecMat.o SAT_Time.o SAT_RefSys.o SAT_Force.o SAT_DE.o Polygon_Relation.o nrlmsise-00_data.o nrlmsise-00.o eopspw.o APC_Planets.o
	g++ -o plan test.o SAT_VecMat.o SAT_Time.o SAT_RefSys.o SAT_Force.o SAT_DE.o Polygon_Relation.o nrlmsise-00_data.o nrlmsise-00.o eopspw.o APC_Planets.o
test.o: test.cpp Polygon_Relation.h SAT_Const.h SAT_DE.h SAT_Force.h SAT_RefSys.h SAT_Time.h SAT_VecMat.h APC_Planets.h eopspw.h GNU_iomanip.h
	g++ -c test.cpp -std=c++11
SAT_VecMat.o: SAT_VecMat.cpp SAT_VecMat.h GNU_iomanip.h
	g++ -c SAT_VecMat.cpp -std=c++11 
SAT_Time.o: SAT_Time.cpp SAT_Time.h GNU_iomanip.h
	g++ -c SAT_Time.cpp -std=c++11 
SAT_RefSys.o: SAT_RefSys.cpp SAT_VecMat.h SAT_RefSys.h SAT_Const.h
	g++ -c SAT_RefSys.cpp -std=c++11 
SAT_Force.o: SAT_Force.cpp SAT_Force.h SAT_Const.h SAT_Time.h SAT_RefSys.h SAT_VecMat.h nrlmsise-00.h eopspw.h
	g++ -c SAT_Force.cpp -std=c++11
SAT_DE.o: SAT_DE.cpp SAT_DE.h
	g++ -c SAT_DE.cpp -std=c++11 
Polygon_Relation.o: Polygon_Relation.cpp Polygon_Relation.h Graph_Var.h GNU_iomanip.h
	g++ -c Polygon_Relation.cpp -std=c++11 
nrlmsise-00_data.o: nrlmsise-00_data.cpp
	g++ -c nrlmsise-00_data.cpp  -std=c++11
nrlmsise-00.o: nrlmsise-00.cpp nrlmsise-00.h
	g++ -c nrlmsise-00.cpp  -std=c++11
eopspw.o: eopspw.cpp eopspw.h SAT_Const.h
	g++ -c eopspw.cpp -std=c++11
APC_Planets.o: APC_Planets.cpp APC_Planets.h SAT_Const.h SAT_VecMat.h
	g++ -c APC_Planets.cpp -std=c++11  
clean:
	rm plan test.o SAT_VecMat.o SAT_Time.o SAT_RefSys.o SAT_Force.o SAT_DE.o Polygon_Relation.o nrlmsise-00_data.o nrlmsise-00.o eopspw.o APC_Planets.o
