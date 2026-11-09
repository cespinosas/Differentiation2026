#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "patrones.h"
#include "basics.h" 

using namespace std;

int main(int argc, char **argv) {
    int seed = atoi(argv[1]); 
    int numgenes = atoi(argv[2]);
    int num_muestreo = atoi(argv[3]);
    int combinaciones = 100;

    Patrones pat;
    Basics bas; 

    int* act_por_red = new int[num_muestreo];
    int* inh_por_red = new int[num_muestreo];
    int* p1_por_red = new int[num_muestreo];
    int* p2_por_red = new int[num_muestreo];
    int* act_inh_por_red = new int[num_muestreo];
    int* reg_mutuas_por_red = new int[num_muestreo];
    int* auto_act_por_red = new int[num_muestreo];
    int* auto_inh_por_red = new int[num_muestreo];
    int* tricpos_por_red = new int[num_muestreo];
    int* tricneg_por_red = new int[num_muestreo];
    int* tasa1_por_red = new int[num_muestreo];
    int* tasa2_por_red = new int[num_muestreo];
    int* tasa3_por_red = new int[num_muestreo];
    int* tasa4_por_red = new int[num_muestreo];
    int* ptokio_por_red = new int[num_muestreo];
    int* ntokio_por_red = new int[num_muestreo];
    int* cj_buenos_por_red = new int[num_muestreo];
    int* cj_malos_por_red = new int[num_muestreo];
    string surtido;
    pat.nacimiento(numgenes);
    pat.nacimientobool(numgenes);

    ofstream act_inh_random("Carpeta_red/act_inh_random.txt");
    act_inh_random << "Activaciones totales" << "\t" << "Inhibiciones totales" << "\t" << "Activaciones mutuas" << "\t" << "Inhibiciones mutuas" << "\t" << "Auto-activaciones" << "\t"
               << "Auto-inhibiciones" << "\t" << "Auto-regulaciones totales" << "\t" << "Activaciones-inhibiciones" << "\t" << "Regulaciones mutuas" << "\n";
    ofstream feedbackl_random("Carpeta_red/3ciclos_random.txt");
    feedbackl_random << "3ciclos positivos" << "\t" << "3ciclos negativos" << "\t" << "Tipo de asa 1" << "\t" << "Tipo de asa 2" << "\t" << "Tipo de asa 3" << "\t" << "Tipo de asa 4" << "\n";
    ofstream total_random("Carpeta_red/total_random.txt");
    total_random << "Positivos en Tokio" << "\t" << "Negativos en Tokio" << "\n";
    ofstream conjugados_random("Carpeta_red/conjugados_random.txt");
    conjugados_random << "Patrones conjugados buenos" << "\t" << "Patrones conjugados malos" << "\n";

    for (int z = 0; z < num_muestreo; z++) {
        pat.limpiar_conjugados(numgenes);
        surtido = "Carpeta_red/Aleatorias/network" + bas.inttostring(z + 1) + ".txt";
        pat.aduana(surtido, numgenes);

        int act = pat.activacionestotales(numgenes);
        int inh = pat.inhibicionestotales(numgenes);
        int p1 = pat.activacionesmutuas(numgenes);
        int p2 = pat.inhibicionesmutuas(numgenes);
        int act_inh = pat.actinh(numgenes);
        int auto_act = pat.autoactivaciones(numgenes);
        int auto_inh = pat.autoinhibiciones(numgenes);
        int auto_tot = pat.totalauto(numgenes);
        int reg_mut = pat.totalmutuas(numgenes);
        pat.triciclos(numgenes, combinaciones);
        pat.conjugados(numgenes);

        act_por_red[z] = act;
        inh_por_red[z] = inh;
        p1_por_red[z] = p1;
        p2_por_red[z] = p2;
        act_inh_por_red[z] = act_inh;
        auto_act_por_red[z] = auto_act;
        auto_inh_por_red[z] = auto_inh;
        reg_mutuas_por_red[z] = reg_mut;
        tricpos_por_red[z] = pat.totaltriciclospositivos();
        tricneg_por_red[z] = pat.totaltriciclosnegativos();
        tasa1_por_red[z] = pat.tipoasa1();
        tasa2_por_red[z] = pat.tipoasa2();
        tasa3_por_red[z] = pat.tipoasa3();
        tasa4_por_red[z] = pat.tipoasa4();
        ptokio_por_red[z] = pat.totalciclospositivos(numgenes);
        ntokio_por_red[z] = pat.totalciclosnegativos(numgenes);
        cj_buenos_por_red[z] = pat.conjugadosbuenos();
        cj_malos_por_red[z] = pat.conjugadosmalos();

        cout << "me lleva la..." << "\n";

        act_inh_random << act << "\t" << inh << "\t" << p1 << "\t" << p2 << "\t" 
                   << auto_act << "\t" << auto_inh << "\t" << auto_tot << "\t" 
                   << act_inh << "\t" << reg_mut << "\n";
                   
        feedbackl_random << tricpos_por_red[z] << "\t" << tricneg_por_red[z] << "\t" 
                   << tasa1_por_red[z] << "\t" << tasa2_por_red[z] << "\t" << tasa3_por_red[z] << "\t" << tasa4_por_red[z] << "\n";

        total_random << ptokio_por_red[z] << "\t" << ntokio_por_red[z] << "\n";
                   
        conjugados_random << cj_buenos_por_red[z] << "\t" << cj_malos_por_red[z] << "\n";
    }

    act_inh_random.close();
    feedbackl_random.close();
    total_random.close();
    conjugados_random.close();

    cout << "en verdad aquí?" << "\n";

    ofstream avg_act_inh("Carpeta_red/act_inh_avgrandom.txt");
    avg_act_inh << "Nodos" << "\t" << "Interacciones" << "\t" << "Activaciones" << "\t" << "Inhibiciones" << "\t" << "Activaciones mutuas (promedio)" << "\t" << "Error estándar" << "\t"
               << "Inhibiciones mutuas (promedio)" << "\t" << "Error estándar" << "\t" << "Regulaciones mutuas (mismo signo)(promedio)" << "\t" << "Error estándar" << "\t"
               << "Activaciones-inhibiciones (promedio)" << "\t" << "Error estándar" << "\t" <<"Auto-activaciones (promedio)" << "\t" << "Error estándar" << "\t"
               << "Auto-inhibiciones (promedio)" << "\t" << "Error estándar" << "\n";
    
    double num_interacciones_avg = bas.get_mean(act_por_red, num_muestreo) + bas.get_mean(inh_por_red, num_muestreo);
    
    avg_act_inh << numgenes << "\t" << num_interacciones_avg << "\t" 
               << bas.get_mean(act_por_red, num_muestreo) << "\t" << bas.get_mean(inh_por_red, num_muestreo) << "\t" //si esto no sale, se finiquita, no lo puedo creer
               << bas.get_mean(p1_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(p1_por_red, num_muestreo) << "\t"
               << bas.get_mean(p2_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(p2_por_red, num_muestreo) << "\t"
               << bas.get_mean(reg_mutuas_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(reg_mutuas_por_red, num_muestreo) << "\t"
               << bas.get_mean(act_inh_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(act_inh_por_red, num_muestreo) << "\t"
               << bas.get_mean(auto_act_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(auto_act_por_red, num_muestreo) << "\t"
               << bas.get_mean(auto_inh_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(auto_inh_por_red, num_muestreo) << "\n";
    avg_act_inh.close();

    ofstream avg_feedbackl("Carpeta_red/3ciclos_avgrandom.txt");
    avg_feedbackl << "3ciclos positivos (promedio)" << "\t" << "Error estándar" << "\t" << "3ciclos negativos (promedio)" << "\t" << "Error estándar" << "\t"
               << "Tipo de asa 1 (promedio)" << "\t" << "Error estándar" << "\t" << "Tipo de asa 2 (promedio)" << "\t" << "Error estándar"
               << "\t" << "Tipo de asa 3 (promedio)" << "\t" << "Error estándar" << "\t" << "Tipo de asa 4 (promedio)" << "\n";
    avg_feedbackl << bas.get_mean(tricpos_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tricpos_por_red, num_muestreo) << "\t"
               << bas.get_mean(tricneg_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tricneg_por_red, num_muestreo) << "\t"
               << bas.get_mean(tasa1_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tasa1_por_red, num_muestreo) << "\t"
               << bas.get_mean(tasa2_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tasa2_por_red, num_muestreo) << "\t"
               << bas.get_mean(tasa3_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tasa3_por_red, num_muestreo) << "\t"
               << bas.get_mean(tasa4_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(tasa4_por_red, num_muestreo) << "\n";
    avg_feedbackl.close();

    ofstream avg_total("Carpeta_red/total_avgrandom.txt");
    avg_total << "Positivos en Tokio (promedio)" << "\t" << "Error estándar" << "\t" << "Negativos en Tokio (promedio)" << "\t" << "Error estándar" << "\n";
    avg_total << bas.get_mean(ptokio_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(ptokio_por_red, num_muestreo) << "\t"
             << bas.get_mean(ntokio_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(ntokio_por_red, num_muestreo) << "\n";
    avg_total.close();

    ofstream avg_conj("Carpeta_red/conjugados_avgrandom.txt");
    avg_conj << "Patrones conjugados buenos (promedio)" << "\t" << "Error estándar" << "\t" << "Patrones conjugados malos (promedio)" << "\t" << "Error estándar" << "\n";
    avg_conj << bas.get_mean(cj_buenos_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(cj_buenos_por_red, num_muestreo) << "\t"
             << bas.get_mean(cj_malos_por_red, num_muestreo) << "\t" << bas.get_sample_stderr(cj_malos_por_red, num_muestreo) << "\n";
    avg_conj.close();

    delete[] act_por_red;
    delete[] inh_por_red;
    delete[] p1_por_red;
    delete[] p2_por_red;
    delete[] act_inh_por_red;
    delete[] reg_mutuas_por_red;
    delete[] auto_act_por_red;
    delete[] auto_inh_por_red;
    delete[] tricpos_por_red;
    delete[] tricneg_por_red;
    delete[] tasa1_por_red;
    delete[] tasa2_por_red;
    delete[] tasa3_por_red;
    delete[] tasa4_por_red;
    delete[] ptokio_por_red;
    delete[] ntokio_por_red;
    delete[] cj_buenos_por_red;
    delete[] cj_malos_por_red;

    pat.liberacionbool(numgenes);
    pat.liberacion(numgenes);

    cout << "deja entendu GOATed" << "\n";

    return 0;
}