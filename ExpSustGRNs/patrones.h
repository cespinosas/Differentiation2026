#ifndef PATRONES_H
#define PATRONES_H
#include <string>
#include "basics.h" 

using std::string;

class Patrones {
private:
    Basics bas; 
    int** adjacency_matrix; 
    int tamanio; 
    int triciclos_positivos;
    int triciclos_negativos;
    int tipodeasa1;
    int tipodeasa2;
    int tipodeasa3;
    int tipodeasa4;
    bool** patrones_conjugados; 
    int num_patronesdistintos;
    int conjugado_bueno;
    int conjugado_malo;

public:
    void nacimiento(int numgenes);
    void aduana(string origen, int numgenes);
    void liberacion(int numgenes);
    void nacimientobool(int numgenes);
    void liberacionbool(int numgenes);
    void limpiar_conjugados(int numgenes);
    
    int activaciones(int numgenes);
    int inhibiciones(int numgenes);
    int activacionestotales(int numgenes);
    int inhibicionestotales(int numgenes);
    int activacionesmutuas(int numgenes);
    int inhibicionesmutuas(int numgenes);
    int actinh(int numgenes);
    int autoactivaciones(int numgenes);
    int autoinhibiciones(int numgenes);
    int totalinteracciones(int numgenes);
    int totalmutuas(int numgenes);
    int totalauto(int numgenes);
    void triciclos(int numgenes, int combinaciones);
    void conjugados(int numgenes);
    int totaltriciclospositivos();
    int totaltriciclosnegativos();
    int totalciclospositivos(int numgenes);
    int totalciclosnegativos(int numgenes);
    int tipoasa1();
    int tipoasa2();
    int tipoasa3();
    int tipoasa4();
    int conjugadosbuenos();
    int conjugadosmalos();
};

#endif