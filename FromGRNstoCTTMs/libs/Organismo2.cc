#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "alea.h"
#include "basics.h"
#include "graphi0.h"
#include "Organismo2.h"
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <queue>
#include <algorithm> 
#include <utility>
#include <sstream>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;
using namespace std;

Organismo_2::Organismo_2(){

}

Organismo_2::Organismo_2(Alea& seed){
    star_rng(seed);
}

void Organismo_2::star_rng(Alea& seed){
    ale = seed;
}

void Organismo_2::make_org(GraphI0& grn2, int genes){

	//string Archivo = "atractores.txt";//SOLO SI NECESITAS IMPRIMIR ATT

    eracciones = 0;
	elev = pow(2,genes);
    grn_org.copy(grn2);
    //grn_org.printnw(cout);

    initial_vector(genes);
    contador();
    num_hileras();
    new_find_all_att(genes);
    diff_org.make_nw(fil, true);
    new_fill_diff_cell_nw(genes);
    //diff_org.printnw(cout);
	//get_all_att(Archivo, genes);

	//delete_grn();//ALGUNAS FUNCIONES NECESITAN QUE NO SE BORREN, REVISARLO ANTES DE USAR
	//delete_diff();


}

void Organismo_2::make_multi_org(GraphI0& grn2, int genes, int arch){

	string Archivo = "output_chunks/atractores_" + bas.inttostring(arch) + ".txt";//SOLO SI NECESITAS IMPRIMIR ATT

    eracciones = 0;
	elev = pow(2,genes);
    grn_org.copy(grn2);
    //grn_org.printnw(cout);

    initial_vector(genes);
    contador();
    num_hileras();
    new_find_all_att(genes);
    diff_org.make_nw(fil, true);
    new_fill_diff_cell_nw(genes);
    //diff_org.printnw(cout);
	get_all_att(Archivo, genes);

	delete_grn();//ALGUNAS FUNCIONES NECESITAN QUE NO SE BORREN, REVISARLO ANTES DE USAR
	delete_diff();


}

void Organismo_2::initial_vector(int genes){

	ve = new int[genes];	
	for(int i=0; i<genes; i++)
		ve[i] = -1;

}

void Organismo_2::contador(){
	
	cont = new int [elev];

}

void Organismo_2::num_hileras(){

	phil = new int[elev];

}

void Organismo_2::new_find_all_att(int genes){

	int i,j,k,y,x,z;
	fil = 0;
	cicatt = 0;
	sies = false;
	yaesta = false;
	allatt = new int **[elev];

	for (k = 0; k < elev; k++){

		grn_org.set_as_state(ve);
		grn_org.new_find_an_attractor();
 		hileras = grn_org.attractor_size();
			
		provatt = new int *[hileras];
		for (i=0; i < hileras; i++){
			provatt[i] = new int [genes];
		}
	
	
		for (i = 0; i < hileras; i++){
		for (j = 0; j < genes; j++)
			provatt[i][j] = grn_org.attractor_element(i,j);  
		}

		if (hileras > 1){
		sies = true;
		}
		  
		grn_org.clear_attractor();
	
		for (z = 0; z < fil; z++){
			yaesta = bas.eqmatrix_rot(provatt, hileras, genes, allatt[z], phil[z], genes);
	
			if (yaesta){
				cont[z]++;
				break;
			}
		}
	
		if (!yaesta){
			allatt[fil] = new int*[hileras];
  			for (j = 0; j < hileras; j++)
    			allatt[fil][j] = new int[genes];
			for (y = 0; y < hileras; y++){
				for (x = 0; x < genes; x++){
				allatt[fil][y][x] = provatt[y][x];
				}
			}
		
					
		cont[fil] = 1;
		phil[fil] = hileras;	
		fil++;
		}
	
		if (k < (elev-1)) {
		bas.next_in_vector(ve, genes, -1, 1, 2);
		}

		for (i=0; i < hileras; i++)
			delete [] provatt[i];
		delete [] provatt; 

	
        if (sies){
            cicatt++;
            sies = false;
        }
		
	
	}
	
}


void Organismo_2::new_fill_diff_cell_nw(int genes){ 

	int a,b,i,j,z,x,k;
	int hilera = 0;
	double peso;
	yaesta = false;
	int **tempatt;
	double var = 0; // periodo en atractor temporal
	
	int ***temporal;
	temporal = new int **[fil];
		for(j = 0; j < fil; j++){
		hilera = phil[j];
		temporal[j] = new int*[hilera];
  			for (i = 0; i < hilera; i++){
    			temporal[j][i] = new int[genes];
			}
		}
		
		
	for (z = 0; z < fil; z++){
		for (i = 0; i < phil[z]; i++){
			for (j = 0; j < genes; j++){
			temporal[z][i][j] = 0;
			}
		}
	}
	
	
	for (i=0; i<fil; i++){
		for(j=0; j<phil[i]; j++){
			for (x=0 ; x<genes; x++){
			temporal[i][j][x] = allatt[i][j][x];
			}
		}
	}	
	
	//cout << "Aqui empiezan a hacerse las comparaciones." << endl;	
	for (i = 0; i < fil; i++){
		for(j = 0; j < phil[i]; j++){
		
		peso = 1/(double)phil[i];
		
		//cout << "--Este es el de prueba que no sirve para nada: --\n";
		grn_org.set_as_state(temporal[i][j]);
		grn_org.new_find_an_attractor();
		//grn_org.print_attractor(cout);
		grn_org.clear_attractor();
		
		
		for (x = 0 ; x < genes; x++){
		temporal[i][j][x] = temporal[i][j][x] * -1;
			//for (e = 0; e < genes; e++)
			//cout << temporal[i][j][e];
			//cout << endl;
			
			//cout << "--Atractor despues de modificar estado: --.\n";
			grn_org.set_as_state(temporal[i][j]);			
			grn_org.new_find_an_attractor();
			//grn.print_attractor(cout);
			var = grn_org.attractor_size();
			
			tempatt = new int*[grn_org.attractor_size()];
			for (k = 0; k < grn_org.attractor_size(); k++){
				tempatt[k] = new int[genes];
			}
			
			for(a = 0; a < grn_org.attractor_size(); a++){
				for(b = 0; b < genes; b++){
				tempatt[a][b] = grn_org.attractor_element(a,b);	
				}
			}
			
			/*for(a = 0; a < grn.attractor_size(); a++){
				for(b = 0; b < genes; b++){
				cout << tempatt[a][b] << " ";	
				}
			cout << "\n";
			}*/
			
			grn_org.clear_attractor();
			
			for (a = 0; a < fil; a++){
			yaesta = bas.eqmatrix_rot(tempatt, var, genes, allatt[a], phil[a], genes);
						
				if(yaesta){
				//cout <<"---Este atractor ya esta.---\n";
				diff_org.chanint(i, a, peso);
				}
				
				if(yaesta && a != i){
				//cout << "---Ademas de que ya esta le sumamos a el cambio pesado.---\n"; 
				//cout << peso << endl;
				//cout << num_wed << endl;
				}
				
			}
			
							
			temporal[i][j][x] = temporal[i][j][x] * -1;
			
			
			}
		
		peso = 0;
	
		}
	}
		
	for (i=0; i < fil; i++){
		for(j=0; j < phil[i]; j++)
			delete [] temporal[i][j];
		delete [] temporal[i];
	}
	delete [] temporal;
	
}

void Organismo_2::delete_grn(){

    grn_org.clear();

}

void Organismo_2::delete_diff(){

    diff_org.clear();

}

void Organismo_2::get_all_att(string Archivo, int genes){

	ofstream fa(Archivo, ios::app);
	if(fa.is_open()){
		
		for (int i = 0; i < fil; i++){
			for (int j = 0; j < phil[i]; j++){
				for (int x = 0; x < genes; x++){
					fa << allatt[i][j][x] << ",";
				}
			}
			fa << "\t";
		}
		fa << endl;

	}
	fa.close();

}

void Organismo_2::make_org_zero(GraphI0& grn2, int genes, int itera, int arch){

	string Archivo2 = "output_chunks/lista_cero" + bas.inttostring(arch) + ".txt";
    eracciones = 0;
	elev = pow(2,genes);
	int val = (arch-1) * 50000;

    grn_org.copy(grn2);
    grn_org.printnw(cout);

    initial_vector(genes);
    contador();
    num_hileras();
    new_find_all_att(genes);
    diff_org.make_nw(fil, true);
    new_fill_diff_cell_nw(genes);
    //diff_org.printnw(cout);
	zero_att_position(Archivo2, genes, val, itera);

	delete_grn();
	delete_diff();


}

void Organismo_2::zero_att_position(string Archivo, int genes, int val, int num_nw){

	ofstream fa(Archivo, ios::app);
	int iguales = 0;
	int comparacion = fil*fil;
  	diff_org_check.make_nw(fil, true);

	
	//------Modificar dependiendo la matriz buscada.--------
	diff_org_check.chanint(0, 0, 1);
	diff_org_check.chanint(1, 0, 1);
	diff_org_check.chanint(0, 1, 0);
	diff_org_check.chanint(1, 1, 0); 

	//diff_org.printnw(cout);
	for(int i = 0; i < fil; i++){
		for(int j = 0; j < fil; j++){
			if(diff_org.nw[j][i] > 0){
				diff_org.nw[j][i] = 1;
			}
		}
	}
	diff_org.printnw(cout);
	diff_org_check.printnw(cout);

	for(int i = 0; i < fil; i++){
		for(int j = 0; j < fil; j++){
			if(diff_org_check.nw[j][i] == diff_org.nw[j][i])
			iguales++;
		}
	}

	//cout << iguales << endl;

	if(iguales == comparacion){
		fa << num_nw + val << endl;
	}
	fa.close();

}