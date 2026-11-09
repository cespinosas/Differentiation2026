#include "patrones.h"
#include <iostream>
#include <fstream>

using namespace std; 

void Patrones::nacimiento(int numgenes) {
    tamanio = numgenes; 
    
    adjacency_matrix = new int*[tamanio];
    for (int j = 0; j < tamanio; j++) {
        adjacency_matrix[j] = new int[tamanio];
    }
    bas.fillmat0(adjacency_matrix, tamanio, tamanio);
}

void Patrones::nacimientobool(int numgenes) {
    num_patronesdistintos = 9;
    patrones_conjugados = new bool*[numgenes];
    for (int i = 0; i < numgenes; i++) {
        patrones_conjugados[i] = new bool[num_patronesdistintos];
    }
    bas.fillmat0(patrones_conjugados, numgenes, num_patronesdistintos);
}

void Patrones::limpiar_conjugados(int numgenes) {
    bas.fillmat0(patrones_conjugados, numgenes, num_patronesdistintos);
}

void Patrones::aduana(string origen, int numgenes) {
    ifstream adj_file(origen.c_str()); 
    if (!adj_file.is_open()) {  
        cerr << "Error al abrir el archivo: " << origen << endl; 
        return;    
    }    

    for (int l = 0; l < numgenes; l++) {    
        for (int n = 0; n < numgenes; n++) {
            adj_file >> adjacency_matrix[l][n];
        }
    }
    adj_file.close();
}

void Patrones::liberacion(int numgenes) {
    if (adjacency_matrix != nullptr) {
        for (int i = 0; i < numgenes; i++) {
            delete[] adjacency_matrix[i];
        }
        delete[] adjacency_matrix;
    }
}

void Patrones::liberacionbool(int numgenes) {
    if (patrones_conjugados != nullptr) {
        for (int i = 0; i < numgenes; i++) {
            delete[] patrones_conjugados[i];
        }
        delete[] patrones_conjugados;
    }
}

int Patrones::activaciones(int numgenes){
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        for (int k = 0; k < numgenes; k++) {
            if (i != k && adjacency_matrix[k][i] > 0) cuenta++;
        }
    }
    return cuenta;
}

int Patrones::inhibiciones(int numgenes){
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        for (int k = 0; k < numgenes; k++) {
            if (i != k && adjacency_matrix[k][i] < 0) cuenta++;
        }
    }
    return cuenta;
}

int Patrones::activacionestotales(int numgenes) {
    return activaciones(numgenes) + autoactivaciones(numgenes);
}

int Patrones::inhibicionestotales(int numgenes) {
    return inhibiciones(numgenes) + autoinhibiciones(numgenes);
}

int Patrones::activacionesmutuas(int numgenes) {
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        for (int k = i + 1; k < numgenes; k++) {
            if (adjacency_matrix[k][i] > 0 && adjacency_matrix[i][k] > 0) {
                cuenta++;
                patrones_conjugados[i][0] = true; 
                patrones_conjugados[k][0] = true;
            }
        }
    }
    return cuenta;
}

int Patrones::inhibicionesmutuas(int numgenes) {
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        for (int k = i + 1; k < numgenes; k++) {
            if (adjacency_matrix[k][i] < 0 && adjacency_matrix[i][k] < 0) {
                cuenta++;
                patrones_conjugados[i][1] = true; 
                patrones_conjugados[k][1] = true;
            }
        }
    }
    return cuenta;
}

int Patrones::actinh(int numgenes) {
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        for (int k = i + 1; k < numgenes; k++) {
            if (adjacency_matrix[k][i] * adjacency_matrix[i][k] < 0) {
                cuenta++; 
                patrones_conjugados[i][2] = true; 
                patrones_conjugados[k][2] = true;
            }
        }
    }
    return cuenta;
}

int Patrones::autoactivaciones(int numgenes) {
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        if (adjacency_matrix[i][i] > 0) {
            cuenta++;
            patrones_conjugados[i][3] = true; 
        }
    }
    return cuenta;
}

int Patrones::autoinhibiciones(int numgenes) {
    int cuenta = 0;
    for (int i = 0; i < numgenes; i++) {
        if (adjacency_matrix[i][i] < 0) {
            cuenta++;
            patrones_conjugados[i][4] = true; 
        }
    }
    return cuenta;
}

int Patrones::totalinteracciones(int numgenes) {
    return activacionestotales(numgenes) + inhibicionestotales(numgenes);
}

int Patrones::totalmutuas(int numgenes) {
    return activacionesmutuas(numgenes) + inhibicionesmutuas(numgenes);
}

int Patrones::totalauto(int numgenes) {
    return autoactivaciones(numgenes) + autoinhibiciones(numgenes);
}

void Patrones::triciclos(int numgenes, int combinaciones) {
    triciclos_positivos = 0;
    triciclos_negativos = 0;
    tipodeasa1 = 0;
    tipodeasa2 = 0;
    tipodeasa3 = 0;
    tipodeasa4 = 0;

    int l, i, k, n, j, z, m, w, q, g;
    int triada = 3;
    int triadas[triada] = {0, 0, 0};
    int signos[triada] = {0, 0, 0};
    int negativos = 0;
    int triadas_ordenadas[triada] = {0, 0, 0};
    int les_gana = 0;
    int positivos = 0;
    int **triadas_juntas;
    bool valio, valio2;
    
    triadas_juntas = new int*[combinaciones];
    for (j = 0; j < combinaciones; j++){
        triadas_juntas[j] = new int[triada];
    }
    bas.fillmat0(triadas_juntas, combinaciones, triada);
    
    int suma_macabrona = 0;
    bool rec_detecta = false;
    bool detecta = false;

    for (i = 0; i < numgenes; i++){
        for (k = 0; k < numgenes; k++){
            if ((adjacency_matrix[i][k] != 0) && (i != k)){
                //cout << "simon" << "\n";
                signos[0] = adjacency_matrix[i][k];
                triadas[0] = i;
                //cout << "triadas[0] es: " << triadas[0] << "\n";
                triadas[1] = k;
                //cout << "triadas[1] es: " << triadas[1] << "\n";

                for (l = 0; l < numgenes; l++){
                    //cout << "simona" << "\n";
                    if ((adjacency_matrix[triadas[1]][l]) != 0 && (adjacency_matrix[l][triadas[0]] != 0) && ((l != triadas[0]) && (l != triadas[1]))){
                        signos[1] = adjacency_matrix[triadas[1]][l];
                        signos[2] = adjacency_matrix[l][triadas[0]];
                        //cout << "los signos: " << signos[0] << "\t" << signos[1] << "\t" << signos[2] << "\n";
                        triadas[2] = l;

                        //cout << "triadas[0] es: " << triadas[0] << "\n";
                        //cout << "triadas[1] es: " << triadas[1] << "\n";
                        //cout << "triadas[2] es: " << triadas[2] << "\n";

                        //if ((triadas[0] == triadas[1]) || (triadas[0] == triadas[2]) || (triadas[1] == triadas[2])){
                          //  cout << "esta ni debería de estar" << "\n";
                        //}
                        //cout << "triadas[0] es: " << triadas[0] << "\n";
                        //cout << "Ahí te voy" << "\n";
                        //for (z=0; z<triada; z++){
                        //cout << triadas[z] << "\t";
                        //}

                        //cout << "worales" << "\n";
                        for (n = 0; n < 3; n++){
                            les_gana = 0;
                            for (z = 0; z < triada; z++){
                                //cout << "epale" << "\n";
                                //cout << "pos" << les_gana << "\n";
                                if (triadas[n] < triadas[z]){
                                    //cout << "gueno" << "\t";
                                    les_gana++;
                                    //cout << les_gana << "\n";
                                } 
                            }
                            triadas_ordenadas[les_gana] = triadas[n];
                            //cout << triadas_ordenadas[les_gana] << "\n";
                        }

                        /*if (i==4 && k ==10 && l ==12){
                            for (z=0; z<triada; z++){
                                cout << triadas_ordenadas[z] << "\t";
                            }
                            cout << "\n" << "\n";
                            for (a=0;a<combinaciones;a++){
                                for (b=0;b<triada;b++){
                                    cout << triadas_juntas[a][b] << " ";
                                }
                                cout << "\n";
                            }
                        } */

                        detecta = false;
                        rec_detecta = false;
                        for (q = 0; q < combinaciones; q++){
                            //suma_valios = 0;
                            //cout << "la suma macabrona va en: " << q << "\n";
                            valio = bas.eqvec(triadas_ordenadas, triada, triadas_juntas[q], triada);
                            if (valio){
                                //suma_valios++;
                                for (g = q + 1; g < combinaciones; g++){
                                    valio2 = bas.eqvec(triadas_ordenadas, triada, triadas_juntas[g], triada);
                                    if (valio2){
                                        //suma_valios++;
                                        rec_detecta = true;
                                        //if (rec_detecta){
                                            //cout << "Pues sí la detecta" << endl;
                                        //}
                                        break;
                                    }
                                }
                                if (rec_detecta) break;
                                if (!rec_detecta && valio){
                                    if (adjacency_matrix[k][i] == 0 || adjacency_matrix[l][k] == 0 || adjacency_matrix[i][l] == 0){
                                        //cout << triadas_juntas[q][0] << "\t" << triadas_juntas[q][1] << "\t" << triadas_juntas[q][2] << " y la q es: " << q << "\n";
                                        //cout << "tons si entra" << "\n";
                                        detecta = true;
                                        //if (detecta){
                                        //cout << "Pues sí la detecta x2" << endl;
                                        //}
                                        //cout << "sí jala" << "\n";
                                        break;
                                    }
                                }
                            }
                        }
                        //cout << "detecta o neh (0 es neh): " << detecta << "\t" << rec_detecta << "\n";
                        if (!detecta && !rec_detecta){
                            positivos = 0;
                            negativos = 0;
                            //cout << "ahi va" << "\n";
                            for (m = 0; m < triada; m++){
                                triadas_juntas[suma_macabrona][m] = triadas_ordenadas[m];
                            }
                            suma_macabrona++;
                            for (m = 0; m < triada; m++){
                                if (signos[m] > 0) positivos++;
                                else negativos++;
                            }
                            if ((negativos%2) == 0) triciclos_positivos++;           
                            else triciclos_negativos++;
                            
                            if (negativos == 0){
                                tipodeasa1++;
                                patrones_conjugados[triadas[0]][5] = true;
                                patrones_conjugados[triadas[1]][5] = true;
                                patrones_conjugados[triadas[2]][5] = true;
                            }
                            if (negativos == 2){
                                tipodeasa2++;
                                patrones_conjugados[triadas[0]][6] = true;
                                patrones_conjugados[triadas[1]][6] = true;
                                patrones_conjugados[triadas[2]][6] = true;
                            }
                            if (negativos == 1){
                                tipodeasa3++;
                                patrones_conjugados[triadas[0]][7] = true;
                                patrones_conjugados[triadas[1]][7] = true;
                                patrones_conjugados[triadas[2]][7] = true;
                            }
                            if (negativos == 3){
                                tipodeasa4++;
                                patrones_conjugados[triadas[0]][8] = true;
                                patrones_conjugados[triadas[1]][8] = true;
                                patrones_conjugados[triadas[2]][8] = true;
                            }
                        }       
                    }
                }
            }
        }
    }

    for (j = 0; j < combinaciones; j++){
        delete[] triadas_juntas[j];
    }
    delete[] triadas_juntas;
}

void Patrones::conjugados(int numgenes) {
    conjugado_bueno = 0;
    conjugado_malo = 0;
    
    for (int i = 0; i < numgenes; i++){
        if ((patrones_conjugados[i][1] || patrones_conjugados[i][6]) && 
            (patrones_conjugados[i][0] || patrones_conjugados[i][3] || patrones_conjugados[i][5])) {
            cout << "Este gen es bueno: " << i << endl;
            conjugado_bueno++;
        }
        if ((patrones_conjugados[i][1] || patrones_conjugados[i][6]) && 
            (patrones_conjugados[i][2] || patrones_conjugados[i][4] || patrones_conjugados[i][7] || patrones_conjugados[i][8])) {
            cout << "Este gen es malo: " << i << endl;
            conjugado_malo++;
        } 
    }
}

int Patrones::totaltriciclospositivos() {
    return triciclos_positivos;
}

int Patrones::totaltriciclosnegativos() {
    return triciclos_negativos;
}

int Patrones::totalciclospositivos(int numgenes) {
    return autoactivaciones(numgenes) + totalmutuas(numgenes) + totaltriciclospositivos();
}

int Patrones::totalciclosnegativos(int numgenes) {
    return autoinhibiciones(numgenes) + actinh(numgenes) + totaltriciclosnegativos();
}

int Patrones::tipoasa1() {
    return tipodeasa1;
}

int Patrones::tipoasa2() {
    return tipodeasa2;
}

int Patrones::tipoasa3() {
    return tipodeasa3;
}

int Patrones::tipoasa4() {
    return tipodeasa4;
}

int Patrones::conjugadosbuenos() {
    return conjugado_bueno;
}

int Patrones::conjugadosmalos() {
    return conjugado_malo;
}
