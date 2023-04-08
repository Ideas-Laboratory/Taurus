g++ -O3 -fopenmp -c src/solver.cpp  src/RepulsiveForce.cpp src/Multilevel.cpp src/layout.cpp  test.cpp -I include/
g++ -O3 -fopenmp Multilevel.o layout.o solver.o RepulsiveForce.o test.o -I include/ -o testc
./testc