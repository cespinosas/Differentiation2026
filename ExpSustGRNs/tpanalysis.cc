#include <iostream>
#include <fstream>
#include <cstdlib>
#include "patrones.h"

using namespace std;

int main(int argc, char **argv) {
    int seed = atoi(argv[1]); 
    int numgenes = atoi(argv[2]);
    int combinaciones = 100; 

    Patrones pat;
    cout << "\n" "\n";

    pat.nacimiento(numgenes);
    pat.nacimientobool(numgenes);

    pat.aduana("Carpeta_red/originalnetwork.txt", numgenes);

    int activaciones = pat.activacionestotales(numgenes);
    int inhibiciones = pat.inhibicionestotales(numgenes);
    int p1 = pat.activacionesmutuas(numgenes);
    int p2 = pat.inhibicionesmutuas(numgenes);
    int act_inh = pat.actinh(numgenes);
    int auto_activaciones = pat.autoactivaciones(numgenes);
    int auto_inhibiciones = pat.autoinhibiciones(numgenes);
    int num_interacciones = pat.totalinteracciones(numgenes);
    int regulaciones_mutuas = pat.totalmutuas(numgenes);
    int auto_totales = pat.totalauto(numgenes);
    pat.triciclos(numgenes, combinaciones);
    pat.conjugados(numgenes);
    int tricpos = pat.totaltriciclospositivos();
    int tricneg = pat.totaltriciclosnegativos();
    int tasa1 = pat.tipoasa1();
    int tasa2 = pat.tipoasa2();
    int tasa3 = pat.tipoasa3();
    int tasa4 = pat.tipoasa4();
    int ptokio = pat.totalciclospositivos(numgenes);
    int ntokio = pat.totalciclosnegativos(numgenes);
    int cj_buenos = pat.conjugadosbuenos();
    int cj_malos = pat.conjugadosmalos();

    cout << "me lleva la..." << "\n";

    cout << "Regulaciones mutuas (mismo signo): " << regulaciones_mutuas << "\n";
    cout << "Activaciones mutuas: " << p1 << "\n"; 
    cout << "Inhibiciones mutuas: " << p2 << "\n";
    cout << "Activaciones-inhibiciones: " << act_inh << "\n";
    cout << "Activaciones totales: " << activaciones << "\n";
    cout << "Inhibiciones totales: " << inhibiciones << "\n";
    cout << "Auto-activaciones totales: " << auto_activaciones << "\n";
    cout << "Auto_inhibiciones totales: " << auto_inhibiciones << "\n";
    cout << "Auto_regulaciones totales: " << auto_totales << "\n";
    cout << "3ciclos positivos: " << tricpos << "\n"; 
    cout << "3ciclos negativos: " << tricneg << "\n";
    cout << "Tipo de asa 1: " << tasa1 << "\n"; 
    cout << "Tipo de asa 2: " << tasa2 << "\n";
    cout << "Tipo de asa 3: " << tasa3 << "\n"; 
    cout << "Tipo de asa 4: " << tasa4 << "\n";
    cout << "Positivos en Tokio: " << ptokio << "\n";
    cout << "Negativos en Tokio: " << ntokio << "\n";
    cout << "Patrones conjugados buenos: " << cj_buenos << "\n"; 
    cout << "Patrones conjugados malos: " << cj_malos << "\n";

    cout << "en verdad aquí?" << "\n";

    ofstream activations_inhibitions("Carpeta_red/act_inh_original.txt");
    activations_inhibitions << "Nodos" << "\t" << "Interacciones" << "\t" << "Activaciones" << "\t" << "Inhibiciones" << "\t" << "Activaciones mutuas" << "\t" << "Inhibiciones mutuas" << "\t"
                            << "Regulaciones mutuas (mismo signo)" << "\t" << "Activaciones-inhibiciones" << "\t" << "Auto-activaciones" << "\t" << "Auto-inhibiciones" << "\n";
    activations_inhibitions << numgenes << "\t" << num_interacciones << "\t" << activaciones << "\t" << inhibiciones << "\t" << p1 << "\t" << p2 << "\t"
                            << regulaciones_mutuas << "\t" << act_inh << "\t" << auto_activaciones << "\t" << auto_inhibiciones << "\n";
    activations_inhibitions.close();

    ofstream feedbackl("Carpeta_red/3ciclos_original.txt");
    feedbackl << "3ciclos positivos" << "\t" << "3ciclos negativos" << "\t" << "Tipo de asa 1" << "\t" << "Tipo de asa 2" << "\t" << "Tipo de asa 3" << "\t" << "Tipo de asa 4" << "\n";
    feedbackl << tricpos << "\t" << tricneg << "\t" << tasa1 << "\t" << tasa2 << "\t" << tasa3 << "\t" << tasa4 << "\n";
    feedbackl.close();

    ofstream total("Carpeta_red/total_original.txt");
    total << "Positivos en Tokio" << "\t" << "Negativos en Tokio" << "\n";
    total << ptokio << "\t" << ntokio << "\n";
    total.close();

    ofstream conjugate_patterns("Carpeta_red/conjugados_original.txt");
    conjugate_patterns << "Patrones conjugados buenos" << "\t" << "Patrones conjugados malos" << "\n";
    conjugate_patterns << cj_buenos << "\t" << cj_malos << "\n";
    conjugate_patterns.close();

    pat.liberacionbool(numgenes);
    pat.liberacion(numgenes);

    cout << "deja entendu GOATed" << "\n";

    return 0;
}