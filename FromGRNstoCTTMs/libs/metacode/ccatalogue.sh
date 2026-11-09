rm makecatalogue makecatalogue.o
g++ -Wall -I/Users/ce/libs/ -c makecatalogue.cc
g++ -o makecatalogue -L/Users/ce/libs/ makecatalogue.o ../basics.o ../alea.o ../graphi.o ../graphc.o ../graphb.o -lgsl -lgslcblas -lm
./makecatalogue 

