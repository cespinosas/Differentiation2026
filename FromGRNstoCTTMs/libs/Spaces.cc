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
#include "Spaces.h"
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

Space::Space(){

}

Space::Space(Alea& seed){

    star_rng(seed);

}

void Space::star_rng(Alea& seed){

    ale = seed;

}

int Space::return_grn(int i, int j){

    /* Te regresa un valor almacenado en grn, usamos esto para que te regrese toda la matriz y almacenarlo 
        en un archivo de texto. Para guarda coherencia que al guardar los valores no esten separados y sea
        una secuencia de números sin espoacios.  */
    
    return grn.obtval(j, i);


}

int Space::return_matriz_nodo(){

    return num_interactions_newnw;

}

int Space::return_nw_positive_interactions(){

    return num_positive_interactions_newnw;

}

void Space::get_n_nw(string Archivo, int filas){

	ifstream nw_file(Archivo.c_str());
	vector<int> num;
	int entry;
	int entry_value = 0;

	while(nw_file >> entry){
		num.push_back(entry);
	}

	grn.make_nw(filas, true);
	for(int i = 0; i < filas; i++){
		for(int j = 0; j < filas; j++){
			grn.fill_nw_value(num[entry_value], j, i);
			entry_value++;
		}
	}
	//grn.printnw(cout);
	nw_file.close();

}

void Space::get_positive_interactions(string Archivo, int i){

    string positive_intera_nw_db;
    int fila = 0;

    ifstream arch(Archivo.c_str());
    if(!arch.is_open()){
        cout << "[ERROR]: NO SE PUDO ABRIR EL ARCHIVO" << endl;
    }

    while(getline(arch, positive_intera_nw_db)){
        if(fila == i){
            break;
        }
        fila++;
    }
	//cout << "[" << positive_intera_nw_db << "]" << endl;
	stringstream ss(positive_intera_nw_db);
	ss >> num_positive_interactions_db;

    arch.close();

}

void Space::get_nwgrn_db(string Archivo, int i, int genes){

    string inw_db;

    int fila = 0;
	int idx_str = 0;
	int signo;
	int num;

    ifstream arch(Archivo.c_str());
    if(!arch.is_open()){
        cout << "[ERROR]: NO SE PUDO ABRIR EL ARCHIVO" << endl;
    }

    while(getline(arch, inw_db)){
        if(fila == i){
            break;
        }
        fila++;
    }

	//cout << "[" << inw_db << "]" << endl;
	grn_check.make_nw(genes, true);

	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			signo = 1;

			if(inw_db[idx_str] == '-'){
				signo = -1;
				idx_str++;
			}

			num = inw_db[idx_str] - '0';
			grn_check.fill_nw_value((signo*num), j, i);
			idx_str++;
		}
	}

	//grn_check.printnw(cout);
    arch.close();

}

void Space::get_nwgrn_original(string Archivo, int i, int genes){

    string inw_db;

    int fila = 0;
	int idx_str = 0;
	int signo;
	int num;

    ifstream arch(Archivo.c_str());
    if(!arch.is_open()){
        cout << "[ERROR]: NO SE PUDO ABRIR EL ARCHIVO" << endl;
    }

    while(getline(arch, inw_db)){
        if(fila == i){
            break;
        }
        fila++;
    }

	//cout << "[" << inw_db << "]" << endl;
	grn.make_nw(genes, true);

	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			signo = 1;

			if(inw_db[idx_str] == '-'){
				signo = -1;
				idx_str++;
			}

			num = inw_db[idx_str] - '0';
			grn.fill_nw_value((signo*num), j, i);
			idx_str++;
		}
	}

	grn.printnw(cout);
    arch.close();

}

void Space::save_grn(string File_name, int genes){

	ofstream archivo(File_name, ios::app); 

    for(int i = 0; i < genes; i++){
        for(int j = 0; j < genes; j++){
            archivo << grn.obtval(j, i);
        }
    }
    archivo << "\n";
    archivo.close();
}

void Space::make_all_nw(const string& nombresal, const string& nombre_pos_intera, const string& nombre_nodo_intera, const string& nombrentr, int n, int start, int end){
	
	/*	Funcion principal a partir de la cual llamamos a todas las demas funciones para
		contruir el espacio de matrices. 	
	 */

	int i,j,z;
	int nw_int = 0;
	int nw_posint = 0;
	num_interactions_newnw = 0;
	int db_size = 1;
	int total_filas_comb = bas.factorial(n);
	save_nw = false;
	one_component = true;
	
	genes = n;
	elev = pow(2, n);
	yaparams = true;

    ofstream fa;
    fa.open(nombresal.c_str());
    ofstream fb;
    fb.open(nombre_nodo_intera.c_str());
    ofstream fc;
    fc.open(nombre_pos_intera.c_str());

	get_n_nw(nombrentr, n);
	//grn.printnw(cout);

	func_combination(n);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(grn.obtval(j,i) != 0){
				nw_int++;
			}
		}
	}
	num_interactions_newnw = nw_int;
	nw_int = 0;

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(grn.obtval(j,i) == 1){
				nw_posint++;
			}
		}
	}

	num_positive_interactions_newnw = nw_posint;
	nw_posint = 0;

    for(int i = 0; i < genes; i++){
        for(int j = 0; j < genes; j++){
            fa << return_grn(i, j); 
        }
    }
    fa << endl;
    fb << return_matriz_nodo() << endl;
    fc << return_nw_positive_interactions() << endl;

 	for(z = start; z < start + end; z++){
		//cout << "Estamos probando la matriz: " << " " << z << endl; 
        save_nw = false;
		one_component = true;
        
		bas.next_in_matrix(grn.nw, n, n, -1, 1, 1);
		//grn.printnw(cout);
		num_interactions_newnw = 0;
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++){
					if(grn.nw[i][j] != 0){
						num_interactions_newnw++;				
					}
				}
			}
		num_positive_interactions_newnw = 0;
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++){
					if(grn.nw[i][j] > 0){
						num_positive_interactions_newnw++;				
					}
				}
			}
			//cout << "Interacciones positivas : " << num_positive_interactions_newnw << endl;

 		check_interaction_rows(n);
    	if(check_decrease_rows){
			check_nodes_interaction(n);
			if(no_empty_rows){
				check_one_component();
				if(one_component){
				for(i = 0; i < db_size; i++){
					get_positive_interactions(nombre_pos_intera, i);
					if(num_positive_interactions_newnw == num_positive_interactions_db){
						get_nwgrn_db(nombresal, i, genes);
						nw_combination(n);
						grn_check.clear();
						if(save_nw == true){
							break;
						}
						num_interactions_bdnw = 0;
					}
				}
				if(save_nw == false){
					for(int i = 0; i < genes; i++){
                        for(int j = 0; j < genes; j++){
                            fa << return_grn(i, j); 
                        }
                    }
					fa << endl;
                    fb << return_matriz_nodo() << endl;
                    fc << return_nw_positive_interactions() << endl;			
					db_size++;
					num_interactions_bdnw = 0;
				}
			    }
		    }   
		}
	} 

	for (i = 0; i < total_filas_comb; i++){
		delete [] nw_combinations[i];
	}
	delete [] nw_combinations;  

    fa.close();
    fb.close();
    fc.close();
			
}

void Space::check_one_component(){

	if(grn.number_of_components() > 1)
		one_component = false;

	grn.clear_comp();

}

void Space::check_interaction_rows(int fil_col){

		int i,j;
		check_decrease_rows = true;
		
		int* num_interaction_rows;
		num_interaction_rows = new int[fil_col];
		for(i = 0; i < fil_col; i++){
			num_interaction_rows[i] = 0;
		}
		
		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				if(grn.nw[i][j] != 0){
					num_interaction_rows[i]++;
				}
			}
		}
		
		for(i = 0; i < fil_col-1; i++){
			if(num_interaction_rows[i] < num_interaction_rows[i+1]){
				check_decrease_rows = false;
			}
		}
		
		delete[] num_interaction_rows;	
	
}

void Space::check_nodes_interaction(int fil_col){

		int i,j;
		no_empty_rows = true;
		
		int* num_interaction_rows;
		num_interaction_rows = new int[fil_col];
		for(i = 0; i < fil_col; i++){
			num_interaction_rows[i] = 0;
		}
		int* num_interaction_col;
		num_interaction_col = new int[fil_col];
		for(i = 0; i < fil_col; i++){
			num_interaction_col[i] = 0;
		}
		
		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				if(grn.nw[i][j] != 0){
					num_interaction_rows[i]++;
				}
			}
		}

		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				if(grn.nw[j][i] != 0){
					num_interaction_col[i]++;
				}
			}
		}
		
		for(i = 0; i < fil_col; i++){
			if(num_interaction_rows[i] == 0 && num_interaction_col[i] == 0){
				no_empty_rows = false;
			}
		}
		
		delete[] num_interaction_rows;	
		delete[] num_interaction_col;
	
}

void Space::nw_combination(int fil_col){
	
	int i,j,k;
	int num_max_col = 0;
    int num_max_entry = 0;
    int nw_review = 0;
    
	num_max_col = bas.factorial(fil_col);
    num_max_entry = pow(fil_col, 2);

/*   	for(int i = 0; i < num_max_col; i++){
		for(int j = 0; j < fil_col; j++){
			cout << nw_combinations[i][j] << " ";
		}
		cout << endl;
	}    */

    int** new_nw;
    new_nw = new int*[fil_col];
    for(i = 0; i < fil_col; i++){
        new_nw[i] = new int[fil_col];
    }
    bas.fillmat0(new_nw, fil_col, fil_col);

	for(k = 0; k < num_max_col; k++){
		nw_review = 0;
		
		for(i = 0; i < fil_col; i++){ 
			for(j = 0; j < fil_col; j++){
			new_nw[nw_combinations[k][i]][nw_combinations[k][j]] = grn.obtval(j,i);
			}
		}
	
		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
			if(new_nw[i][j] == grn_check.obtval(j,i))
			nw_review++;
			}
		}
		
		if(nw_review == num_max_entry){
			save_nw = true;
			break;
		}
	}

	for (i = 0; i < fil_col; i++){
		delete [] new_nw[i];
	}
	delete [] new_nw;  
  
}

void Space::func_combination(int fil_col){

	int num_max_col = 0;
	vector<int> combinacion;
	vector<bool> usadas(fil_col, false);
	int profundidad = 0;
	int fila_actual = 0;
    
	num_max_col = bas.factorial(fil_col);
    
	vector<int> vec_numeros;
	for(int i = 0; i < fil_col; i++){
		vec_numeros.push_back(i);
	}
	vector<int> indices;
	for(int i  = 0; i < fil_col; i++){
		indices.push_back(0);
	}
	
	//int** nw_combinations;
	nw_combinations = new int*[num_max_col];
	for(int i = 0; i < num_max_col; i++){	
		nw_combinations[i] = new int[fil_col];
	}
	bas.fillmat0(nw_combinations, num_max_col, fil_col);

	while (true) {
        if (profundidad == fil_col){
            for (int i = 0; i < fil_col; i++){
                nw_combinations[fila_actual][i] = combinacion[i];
            }
			fila_actual++;

            profundidad--;
            if (profundidad < 0) break;
            usadas[combinacion.back()] = false;
            combinacion.pop_back();
            indices[profundidad]++;
        } else if (indices[profundidad] >= fil_col) {
            indices[profundidad] = 0;
            if (profundidad == 0) break;
            profundidad--;
            usadas[combinacion.back()] = false;
            combinacion.pop_back();
            indices[profundidad]++;
        } else {
            int candidato = vec_numeros[indices[profundidad]];
            if (!usadas[candidato]) {
                combinacion.push_back(candidato);
                usadas[candidato] = true;
                profundidad++;
            } else {
                indices[profundidad]++;
            }
        }
    }

/*    	for(int i = 0; i < num_max_col; i++){
		for(int j = 0; j < fil_col; j++){
			cout << nw_combinations[i][j] << " ";
		}
		cout << endl;
	}   */


 }

void Space::check_nwdb_interactions(const string& nombreDB, int fil_col, int id){

	int i,j;
	num_interactions_bdnw = 0;
	
	get_nwgrn_db(nombreDB, id, genes);

		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				if(grn_check.nw[i][j] != 0){
					num_interactions_bdnw++;
				}
			}
		}
		
/* 		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				cout << grn_check.nw[i][j] << " ";
			}
			cout << endl;
		}
		cout  << endl;
		
	cout << " Numero de interacciones:" <<  " " << num_interactions_bdnw << endl; 
 */		
}

void Space::review_tar_nw(string nombreDB, int db_size, int genes){
	
	int num_feedback = 0;
	string Archivo = "Lista_feedback.txt";
	ofstream fa;
    fa.open(Archivo.c_str());

	for(int i = 0; i < db_size; i++){

		cout << "Numero de red:" << i << endl;
		num_feedback = 0;
		get_nwgrn_original(nombreDB, i, genes);
		grn.printnw(cout);

		for(int j = 0; j < genes; j++){
			for(int k = 0; k < genes; k++){
				if(j == k && grn.nw[k][j] != 0){
				for(int z = 0; z < genes; z++){
					if(grn.nw[k][z] == 1 && grn.nw[z][k] == 1 && k != z){
					num_feedback++;
					}
				}
				}
			}
		}
		cout << num_feedback << endl;
		fa << num_feedback << endl;
	}


}

void Space::get_distance_db(string nombreDB, int db_size, int genes, int att){

	// ------  Define el objetivo, necesario cambiar en cada ejemplo ------ //
	
	string itera;
 	int** tmp_nw;
	tmp_nw = new int*[att];
	for(int i = 0; i < att; i++){	
		tmp_nw[i] = new int[att];
	}
	bas.fillmat0(tmp_nw, att, att);

	diff_check.make_nw(att, true);
	bas.fillmat0(diff_check.nw, att, att);
    diff_check.fill_nw_value(1, 0, 0);
    diff_check.fill_nw_value(1, 1, 1);
    //diff_check.fill_nw_value(1, 2, 2);
    diff_check.fill_nw_value(1, 0, 1);
    //diff_check.fill_nw_value(1, 0, 2); 

	func_combination(att);
	numero_targets = 0; 
	int numero_3 = 0;

	for(int i  = 0;  i < db_size; i++){
		//cout << "ESTAMOS PROBANDO LA RED" << i << endl;
		get_nwgrn_original(nombreDB, i, genes);
		//grn.printnw(cout);
		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
 		if(org2.fil == 2){
			numero_3++;
			//cout << "Numero con 3" << numero_3 << endl;
			//diff.printnw(cout);
			for(int i = 0; i < att; i++){
				for(int j = 0; j < att; j++){
					if(diff.obtval(j, i) != 0)
					tmp_nw[i][j] = 1;
				}
			}
			
			nw_combination_dif(tmp_nw, att);
			if(nw_target){
				//itera = to_string(i + 1);
				bas.inttostring(i+1);
				grn.print_dot_edgecolor("Figuras/regulacion" + itera);
				bas.get_dot_fig("Figuras/regulacion" + itera);
				cout <<  /* "El target es la red"  << */ i + 1 << endl;
			}
			bas.fillmat0(tmp_nw, att, att);
		}
		org2.delete_grn();
		org2.delete_diff();

	}

	//cout << "Numero de targets" << numero_targets << endl;
	//cout << "Numero de redes mayor o igual a 3: " << numero_3 << endl;
	
	int borrar = bas.factorial(att);
	for (int i = 0; i < borrar; i++){
		delete [] nw_combinations[i];
	}
	delete [] nw_combinations;  

}

void Space::comp_all_space(string nombreDB, int genes){
	
	int space_size = 19683;
	//int space_size = 100000;
	string itera;
	string Archivo = "Guardados.txt";

	/* grn.make_nw(genes, true);
	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			grn.fill_nw_value(-1, i, j);
		}
	} */
	get_n_nw(nombreDB, genes);
	grn.printnw(cout);
	
	int att = 2;
	int** tmp_nw;
	tmp_nw = new int*[att];
	for(int i = 0; i < att; i++){	
		tmp_nw[i] = new int[att];
	}
	bas.fillmat0(tmp_nw, att, att);

	diff_check.make_nw(att, true);
	bas.fillmat0(diff_check.nw, att, att);
	//diff_check.fill_nw_value(1, 2, 2); 
    diff_check.fill_nw_value(1, 1, 1);
    //diff_check.fill_nw_value(1, 0, 0);
    diff_check.fill_nw_value(1, 0, 1);
    //diff_check.fill_nw_value(1, 0, 2);

	func_combination(att);
	numero_targets = 0; 

	for(int i = 1; i < space_size+1; i++){
		cout << "Numero de iteracion: " << i << endl;
		bas.next_in_matrix(grn.nw, genes, genes, -1, 1, 1);
		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		if(org2.fil == att){
		//diff.printnw(cout);	
		for(int i = 0; i < att; i++){
			for(int j = 0; j < att; j++){
				if(diff.obtval(j, i) != 0)
				tmp_nw[i][j] = 1;
			}
		}
		nw_combination_dif(tmp_nw, att);
		if(nw_target){
		save_grn(Archivo, genes);
		}
		bas.fillmat0(tmp_nw, att, att);

		}
		org2.delete_grn();
		org2.delete_diff();

	}


}

void Space::ofunc_comb(int fil){//NO FUNCIONA PERO COMO TAMPOCO SE USA CORREGIDA EN SIM_A

/*     int num_max_col = (fil) * (fil-1) * (fil-2);
    int fila_actual = 0;
    int contador = 0;

    vector<int> combinacion;
    vector<bool> usadas1(fil, false);

	nw_combinations = new int*[num_max_col];
	for(int i = 0; i < num_max_col; i++){	
		nw_combinations[i] = new int[fil];
	}
	bas.fillmat0(nw_combinations, num_max_col, fil); */

/*     for(int i = 0; i < num_max_col; i++){
        for(int j = 0; j < fil; j++){
            cout << nw_combinations[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl; */

 /*    for(int i = 0; i < fil; i++){
        if(usadas1[i] == false){
            usadas1[i] = true;
            combinacion.push_back(i);

            for(int j = 0; j < fil; j++){
                if(usadas1[j] == false){
                    usadas1[j] = true;
                    combinacion.push_back(j);

                    for(int z = 0; z < fil; z++){
                        if(usadas1[z] == false){
                            contador++;
                            usadas1[z] = true;
                            combinacion.push_back(z);
                            for(int k = 0; k < combinacion.size(); k++){
                               nw_combinations[fila_actual][combinacion[k]] = k;
                            }
                            combinacion.pop_back();
                            fila_actual++;
                            usadas1[z] = false;
                        }
                    }
                    combinacion.pop_back();
                    usadas1[j] = false;
                }
            }
            combinacion.pop_back();
            usadas1[i] = false;
        }
    } */

/*     for(int i = 0; i < num_max_col; i++){
        for(int j = 0; j < fil; j++){
            cout << nw_combinations[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Total de combinaciones registradas: " << contador << endl; */

}

void Space::nw_combination_dif(int** nw, int fil_col){
	
	int i,j,k;
	int num_max_col = 0;
    int num_max_entry = 0;
    int nw_review = 0;
	nw_target = false;
    
	num_max_col = bas.factorial(fil_col);
    num_max_entry = pow(fil_col, 2); 

    int** new_nw;
    new_nw = new int*[fil_col];
    for(i = 0; i < fil_col; i++){
        new_nw[i] = new int[fil_col];
    }
    bas.fillmat0(new_nw, fil_col, fil_col);

	for(k = 0; k < num_max_col; k++){
		nw_review = 0;
		
		for(i = 0; i < fil_col; i++){ 
			for(j = 0; j < fil_col; j++){
			new_nw[nw_combinations[k][i]][nw_combinations[k][j]] = diff_check.obtval(j,i);
			}
		}
	
		for(i = 0; i < fil_col; i++){
			for(j = 0; j < fil_col; j++){
				if(new_nw[i][j] == nw[i][j])
				nw_review++;
			}
		}

		//cout << "Numero de entradas iguales" << nw_review << endl;

		if(nw_review == num_max_entry){
			nw_target = true;
			//cout << "Esta es el target" << endl;
			numero_targets++;
			break;
		}
	
	}

	for (i = 0; i < fil_col; i++){
		delete [] new_nw[i];
	}
	delete [] new_nw;  
  
}

void Space::comp_two_db(string nombreDB1, string nombreDB2, string nombresal, int size_db1, int size_db2, int fil_col){

	save_nw = false;
    ofstream fa;
    fa.open(nombresal.c_str());
	int total_filas_comb = bas.factorial(fil_col);
	func_combination(fil_col);

	bool borramos_db1 = del_1row(nombreDB1, fil_col);
	bool borramos_db2 = del_1row(nombreDB2, fil_col);

	if(borramos_db1){
		//cout << "BORRAMOS" << endl;
		size_db1 = size_db1-1;
	}
	if(borramos_db2){
		//cout << "BORRAMOS" << endl;
		size_db2 = size_db2-1;
	}
	int size_db3 = size_db1;

	for(int i = 0; i < size_db1; i++){

		get_nwgrn_original(nombreDB1, i, fil_col);
		for(int j = 0; j < fil_col; j++){
			for(int z = 0; z < fil_col; z++){
				fa << return_grn(j, z);
			}
		}
		fa << endl;

	}
	//cout << "Termino de copiar la primera base" << endl;

	for(int i = 0; i < size_db2; i++){
		//cout << "Esta probando la red " << i << " de la segunda base" << endl;
		save_nw = false;
		//get_nwgrn_db(nombreDB2, i, genes);
		get_nwgrn_original(nombreDB2, i, fil_col);
		//grn.printnw(cout);
		for(int j = 0; j < size_db3; j++){
			//cout << "Comparamos con la red " << j << endl;
			//get_nwgrn_original(nombresal, j, fil_col);
			get_nwgrn_db(nombresal, j, fil_col);
			nw_combination(fil_col);
			grn_check.clear();
			if(save_nw == true){
				break;	
			}
		}
		if(save_nw == false){
			//cout << "No se encontro isomorfa por lo que se guarda la red " << i << endl;
			for(int i = 0; i < fil_col; i++){
				for(int j = 0; j < fil_col; j++){
					fa << grn.obtval(i, j); 
				}
			}
			fa << endl;
			size_db3++;
		}
	}

	for (int i = 0; i < total_filas_comb; i++){
		delete [] nw_combinations[i];
	}
	delete [] nw_combinations;  

    fa.close();

}

bool Space::del_1row(string nombreDB, int genes){//ESTA MAL, CORRGIR ORDEN DE BORRAR Y GUARDAR

	string linea;
	bool primera_linea = true;
	bool delete_1row = false;

	ifstream arch;
	arch.open(nombreDB.c_str());
	ofstream temporal;
	temporal.open("temporal.txt");

	one_component = true;

	get_nwgrn_original(nombreDB, 0, genes);
	check_one_component();
	//cout << "componente" << one_component << endl;
	check_interaction_rows(genes);
	//cout << "intera_rows" << check_decrease_rows << endl;

	if(!one_component || !check_decrease_rows){
	cout << "Borramos" << endl;
		while(getline(arch, linea)){
			if(primera_linea){
			primera_linea = false;
			delete_1row = true;
			}
			else{
			temporal << linea << '\n';
			}
		}
		arch.close();
		temporal.close();
		remove(nombreDB.c_str());
		rename("temporal.txt", nombreDB.c_str());
	}else{
	cout << "No_Borramos" << endl;
	arch.close();
	temporal.close();	
	}

	return delete_1row;
	//arch.close();
	//temporal.close();

}

void Space::red_db(string nombreDB, int part_size, int genes, int itera){
	
	int begin = itera*part_size;
	int end = (itera+1)*part_size;
	//cout << begin << end << endl;
	int nw_actual = 0;
	int actual_size = 1;
	int total_filas_comb = bas.factorial(genes);

	string final_name = "salida_combinada_"+bas.inttostring(itera)+".txt";
	func_combination(genes);

    ofstream fa;
    fa.open(final_name.c_str());

	get_nwgrn_original(nombreDB, nw_actual, genes);//si empieza con cero esta bien
	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			fa << return_grn(i, j);
		}
	}
	fa << endl;

	for(int i = begin+1; i < end; i++){
		//cout << "Iteracion_actual:" << i << endl;
		save_nw = false;
		get_nwgrn_original(nombreDB, i, genes);
		for(int j = 0; j < actual_size; j++){
			get_nwgrn_db(final_name, j, genes);
			nw_combination(genes);
			grn_check.clear();
			if(save_nw == true){
			break;
			}
		}
		if(save_nw == false){
		actual_size++;	
		for(int i = 0; i < genes; i++){    
		    for(int j = 0; j < genes; j++){
                fa << return_grn(i, j); 
            }
        }
		fa << endl;
		}
		
		grn.clear();
	}

	for (int i = 0; i < total_filas_comb; i++){
		delete [] nw_combinations[i];
	}
	delete [] nw_combinations;
	
	fa.close();

}

void Space::make_random_wells(){

	int genes = 4;
	int new_value;
	string Archivo = "prueba2.txt";

	grn.make_nw(genes, true);
	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			new_value = ale.randint(-1,2);
			grn.fill_nw_value(new_value, j, i);
		}
	}
	save_grn(Archivo, genes);
    grn.printnw(cout);

	for(int i = 0; i < 235; i++){

		row1 = ale.randint(0, genes);
		row2 = ale.randint(0, genes);
		new_value = ale.randint(-1, 2);
		while(new_value == grn.obtval(row2, row1)){
			new_value = ale.randint(-1, 2);
		}       
		grn.fill_nw_value(new_value, row2, row1);
		save_grn(Archivo, genes);

	}

}

void Space::space_wells_txt(string nombreDB1, string nombreDB2, int db_size, int genes){

	fill_false(nombreDB1, nombreDB2);
	idx_nw = 1;

	for(int i = 1; i < db_size+ 1; i++){
		cout << "Tamaño del pozo:" << " " << well_size << endl;;
		cout << "Numero de pozos:" << " " << i << endl;
		get_nwgrn_original(nombreDB2, idx_nw, genes);
		change_state(nombreDB2, i);
		well_size = 1;

		for(int j = 0; j < 20000; j++){
			//cout << "Mutacion:" << j << endl;
			random_mutation(genes);
			check_mutation(nombreDB2, db_size, genes);
			if(save_nw){
				change_state(nombreDB2, idx_target);
			}
			else if(!save_nw){
				grn.fill_nw_value(saved_value, row2, row1);
			}
			
		}
	}

}

void Space::space_wells(string nombreDB1, string nombreDB2, int db_size, int genes){

	int num_columns = genes*genes;
	int** all_nw;
	all_nw = new int*[db_size];
	for(int i = 0; i < db_size; i++){
		all_nw[i] = new int[num_columns+1];
	}
	bas.fillmat0(all_nw, db_size, num_columns+1);

	int column;
	for(int i = 0; i < db_size; i++){
		get_nwgrn_original(nombreDB1, i, genes);
		column = 0;
		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				all_nw[i][column] = grn.nw[j][z];
				column++;
			}
		}
		grn.clear();
	}

	idx_nw = 1;
	int well_number = 1;
	int temp_well;
	int valor;
	flag = (genes * genes);
	bool aun_redes = true;
	for(int i = 0; i < db_size; i++){
		cout << "Iteracion" << i << endl;
 		valor = -1;
		for(int h = 0; h < db_size; h++){
			if(all_nw[h][flag] == 0){
				valor = h;
				cout << "Matriz numero:" << " " << valor << endl;
				break;
			}
		}
		if(valor == -1){
			break;
		} 
		get_nwgrn_memory(all_nw, valor, genes);
		all_nw[valor][flag] = well_number;

		for(int j = 0; j < 200000; j++){
			//cout << "Mutacion:" << j << endl;
			random_mutation(genes);
 			review_mutation(all_nw, db_size, genes);
			if(save_nw && all_nw[idx_target][flag] == 0 ){
				all_nw[idx_target][flag] = well_number;
			}
 			else if(save_nw && all_nw[idx_target][flag] != 0 && all_nw[idx_target][flag] != well_number){
				temp_well = well_number;
				well_number = all_nw[idx_target][flag];
				cout << "Nuevo pozo:" << well_number << endl;;
				change_all_well(all_nw, temp_well, well_number, db_size);
				well_number = temp_well;
				break;
			} 
			else if(!save_nw){
				grn.fill_nw_value(saved_value, row2, row1);
			} 
			
		}
		well_number++;

	}

	for(int i = 0; i < db_size; i++){
		for(int j = 0; j < num_columns+1; j++){
			cout << all_nw[i][j] << " ";
		}
		cout << endl;
	}

	for (int i = 0; i < genes+1; i++){
		delete [] all_nw[i];
	}
	delete [] all_nw;


}

void Space::random_walks(string nombreDB, int genes, int itera){

	string File_name = "prueba.txt";

	int total_interaction = 9;//Cambiar con quien se compare
	int iteration_interaction = 0;
	grn.make_nw(genes, true);
	get_nwgrn_original(nombreDB, 0, genes);
	cout << "Esta es la red inicial" << endl;
	grn.printnw(cout);

	diff_check.make_nw(3, true);
	diff_check.printnw(cout); 
	//diff_check.change_interaction(0, 0, 1);//Quitar si buscas con o sin oreja
	diff_check.change_interaction(0, 1, 1);
	diff_check.change_interaction(0, 2, 1);
	diff_check.change_interaction(1, 1, 1);
	diff_check.change_interaction(2, 2, 1);
	cout << "Esta es la red target" << endl;
	diff_check.printnw(cout); 

	org2.make_org(grn, genes);
	diff.copy(org2.diff_org);
	diff.printnw(cout);

	for(int x = 0; x < itera; x++){

		cout << "Iteracion:" << x << endl;
		random_mutation(genes);
		grn.printnw(cout);
		iteration_interaction = 0;

		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		cout << "Nueva red de dif:" << endl;
		diff.printnw(cout);

		if(org2.fil == 3 && org2.cicatt == 0){
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] > 0)
					diff.nw[i][j] = 1;
				}
			}
			diff.printnw(cout);

			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] == diff_check.nw[i][j])
					iteration_interaction++;
				}
			}
		}		

		if(total_interaction != iteration_interaction){
			cout << "No son iguales" << endl;
			grn.fill_nw_value(saved_value, row2, row1);
		}
		if(total_interaction == iteration_interaction){
			cout << "Si son iguales" << endl;
			save_grn(File_name, genes);
		}


	}

}

void Space::random_walks_2ears(string nombreDB, int genes, int itera){

	//	Que la funcion no busque una sola, sino los dos patrones de bifurcacion, puede pasar oreja no oreja sin modificar el patron

	string File_name = "prueba.txt";

	int total_interaction = 9;//Cambiar con quien se compare
	int iteration_interaction = 0;
	//int ear_interaction = 0;//Solo si se buscan puentes
	grn.make_nw(genes, true);
	get_nwgrn_original(nombreDB, 0, genes);
	cout << "Esta es la red inicial" << endl;
	grn.printnw(cout);

	diff_check.make_nw(3, true);
	diff_check.printnw(cout); 
	diff_check.change_interaction(0, 0, 1);//Quitar si no se busca puentes
	diff_check.change_interaction(0, 1, 1);
	diff_check.change_interaction(0, 2, 1);
	diff_check.change_interaction(1, 1, 1);
	diff_check.change_interaction(2, 2, 1);
	cout << "Esta es la red target" << endl;
	diff_check.printnw(cout); 

	//--------- Quitar si no se busca puentes------------
	/* diff_ear.make_nw(3, true);
	diff_ear.printnw(cout);
	diff_ear.change_interaction(0, 0, 1);
	diff_ear.change_interaction(0, 1, 1);
	diff_ear.change_interaction(0, 2, 1);
	diff_ear.change_interaction(1, 1, 1);
	diff_ear.change_interaction(2, 2, 1);
	cout << "Esta red tiene oreja" << endl;
	diff_ear.printnw(cout); */

	org2.make_org(grn, genes);
	diff.copy(org2.diff_org);
	diff.printnw(cout);

	for(int x = 0; x < itera; x++){

		cout << "Iteracion:" << x << endl;
		random_mutation(genes);
		grn.printnw(cout);
		iteration_interaction = 0;
		//ear_interaction = 0; //Solo si se buscan puentes

		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		cout << "Nueva red de dif:" << endl;
		diff.printnw(cout);

		if(org2.fil == 3){
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] > 0)
					diff.nw[i][j] = 1;
				}
			}
			diff.printnw(cout);

			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] == diff_check.nw[i][j])
					iteration_interaction++;
				}
			}

			//--------- Quitar si no se busca puentes------------
			/* for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] == diff_ear.nw[i][j])
					ear_interaction++;
				}
			}
			if(ear_interaction == genes*genes){
				cout << "Esta esta a un lado!!!!!!" << endl;
				break;
			} */
		}		

		if(total_interaction != iteration_interaction){
			cout << "No son iguales" << endl;
			grn.fill_nw_value(saved_value, row2, row1);
		}
		if(total_interaction == iteration_interaction){
			cout << "Si son iguales" << endl;
			save_grn(File_name, genes);
		}


	}

}

void Space::random_nw(string nombreDB, int genes, int itera){

	//string File_name = "prueba.txt";
	int valor;

	int total_interaction = 9;//Cambiar con quien se compare
	int iteration_interaction = 0;
	grn.make_nw(genes, true);

	diff_check.make_nw(3, true);
	diff_check.printnw(cout); 
	diff_check.change_interaction(0, 0, 1);
	diff_check.change_interaction(0, 1, 1);
	diff_check.change_interaction(0, 2, 1);
	diff_check.change_interaction(1, 1, 1);
	diff_check.change_interaction(2, 2, 1);
	cout << "Esta es la red target" << endl;
	diff_check.printnw(cout); 

	for(int x = 0; x < itera; x++){

		cout << "Iteracion:" << x << endl;
		for(int i = 0; i < genes; i++){
			for(int j = 0; j < genes; j++){
				valor = ale.randint(-1, 2);
				grn.nw[i][j] = valor;
			}
		}
		grn.printnw(cout);
		iteration_interaction = 0;

		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		cout << "Nueva red de dif:" << endl;
		diff.printnw(cout);

		if(org2.fil == 3 && org2.cicatt == 0){
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] > 0)
					diff.nw[i][j] = 1;
				}
			}
			diff.printnw(cout);

			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if(diff.nw[i][j] == diff_check.nw[i][j])
					iteration_interaction++;
				}
			}
		}		

		if(total_interaction != iteration_interaction){
			cout << "No son iguales" << endl;
			//grn.fill_nw_value(saved_value, row2, row1);
		}
		if(total_interaction == iteration_interaction){
			cout << "Si son iguales" << endl;
			save_grn(nombreDB, genes);
		}


	}

}

void Space::change_all_well(int **mat, int well_number_seek, int new_well_number, int db_size){

	for(int i = 0; i < db_size; i++){
		if(mat[i][flag] == well_number_seek){
			mat[i][flag] = new_well_number;
		}
	}
}

void Space::get_nwgrn_memory(int **mat, int fila, int genes){

	grn.make_nw(genes, true);
	int colum = 0;

	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			grn.nw[i][j] = mat[fila][colum];
			colum++;
		}
	}
	//cout << "GRN" << endl;
	//grn.printnw(cout);

}

void Space::get_nwgrncheck_memory(int **mat, int fila, int genes){

	grn_check.make_nw(genes, true);
	int colum = 0;

	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			grn_check.nw[i][j] = mat[fila][colum];
			colum++;
		}
	}

	//cout << "GRN_CHECK" << endl;
	//grn_check.printnw(cout);

}

void Space::review_mutation(int **mat, int db_size, int genes){

	int sum;
	save_nw = false;
	int total_entry = genes * genes;

	for(int i = 0; i < db_size; i++){
		get_nwgrncheck_memory(mat, i, genes);
		sum = 0;
		for(int k = 0; k < genes; k++){
			for(int j = 0; j < genes; j++){
			if(grn.nw[k][j] == grn_check.nw[k][j])
			sum++;	
			}
		}
		//cout << "sum:" << sum << endl;
		grn_check.clear();

		if(sum == total_entry){
			save_nw = true;
			idx_target = i;
			//label_target = mat[i][flag];
			//cout << "label_target:" << label_target << endl;
			break;
		}

	}

}

void Space::fill_false(string nombreDB, string nombreDB2){

	ifstream input(nombreDB);
    ofstream output(nombreDB2);

    string line;
    while (getline(input, line)) {
        if (!line.empty()) {
            output << line << " 0\n";
        }
    }

    input.close();
    output.close();
}

void Space::change_state(string nombreDB, int target_line){

	fstream arch(nombreDB, ios::in | ios::out);

	int current_line = 1;
	string line;
    while (getline(arch, line)){
        if(current_line == target_line){
            if(!line.empty() && line.back() == '0'){
                streampos position = arch.tellg();
                arch.seekp(position - streamoff(2));
                arch.put('1');
				well_size++;
            }
            break;
        }
        current_line++;
    }

    arch.close();

}

void Space::random_mutation(int genes){ 

	row1 = ale.randint(0, genes);
   row2 = ale.randint(0, genes);
	int new_value;
	new_value = ale.randint(-1, 2);
	//cout << new_value << endl;

    int i = row1;
    int j = row2;
	while(new_value == grn.obtval(j, i)){
		new_value = ale.randint(-1, 2);
		//cout << new_value << endl;
	}
    saved_value = grn.obtval(j, i);       
    grn.fill_nw_value(new_value, j, i);

	//grn.printnw(cout);

}

void Space::check_mutation(string nombreDB, int db_size, int genes){
	
	int sum;
	save_nw = false;
	int total_entry = genes * genes;

	for(int i = 0; i < db_size; i++){
		get_nwgrn_db(nombreDB, i, genes);
		sum = 0;
		for(int k = 0; k < genes; k++){
			for(int j = 0; j < genes; j++){
			if(grn.nw[k][j] == grn_check.nw[k][j])
			sum++;	
			}
		}
		grn_check.clear();
		if(sum == total_entry){
			save_nw = true;
			idx_target = i;
			break;
		}
	}

}

void Space::get_next_idx(string nombreDB, int db_size){

	fstream arch(nombreDB, ios::in | ios::out);

	int current_line = 1;
	string line;
    while (getline(arch, line)){
        if(!line.empty() && line.back() != '0'){
            idx_nw = current_line;
			break;
        }
    current_line++;
	}

    arch.close();

}

void Space::count_one_step(string nombreDB, string nombreDB2, int genes, int db_size1, int db_size2){

	for(int i = 0; i < db_size1; i++){
		get_nwgrn_original(nombreDB, i, genes);
 		for(int j = 0; j < db_size2; j++){
			get_nwgrn_db(nombreDB2, j, genes);
			number_of_diff(genes);
			if(one_step == true){
				cout << i << j << endl;
			}
		} 
	}

}

void Space::number_of_diff(int genes){

	number_diff = 0;
	one_step = false;
	for(int i = 0; i < genes; i++){
		for(int j = 0; j< genes; j++){
			if(grn.nw[i][j] != grn_check.nw[i][j]){
				number_diff++;
			}
		}
	}
	//cout << number_diff << endl;
	if(number_diff == 1){
		one_step = true;
	}

}

void Space::count_number_differences(string nombreDB, string nombreDB2, int genes, int db_size1, int db_size2){
	
	int menor_valor_dif = 999999;//Algo muy grande para que el primero se acepte siempre
	int conteo_menor_valor = 0;
	int suma = 0; 
	for(int i = 0; i < db_size1; i++){
		cout << "Iteracion" << i << endl;
		get_nwgrn_original(nombreDB, i, genes);
 		for(int j = 0; j < db_size2; j++){
			get_nwgrn_db(nombreDB2, j, genes);
			
			suma = 0;
			for(int k = 0; k < genes; k++){
				for(int z = 0; z < genes; z++){
					if(grn.nw[k][z] != grn_check.nw[k][z]){
						suma++;
					}
				}
			}
			if(suma < menor_valor_dif){
				menor_valor_dif = suma;
				conteo_menor_valor = 1;
			}else if (suma == menor_valor_dif){
				conteo_menor_valor++;
			}

			if(suma == 1 && menor_valor_dif == 1){
				cout << "¡¡¡¡¡¡AQUIIIIIII!!!" << conteo_menor_valor << endl;
				grn.printnw(cout);
				grn_check.printnw(cout);
				for(int k = 0; k < genes; k++){
					for(int z = 0; z < genes; z++){
					cout << grn_check.nw[k][z];
					}
				}
				cout << endl;
			}

		} 
	}
	cout << menor_valor_dif << endl;
	cout << conteo_menor_valor << endl;

}

void Space::get_1fb_clusters(string nombreDB, int db_size, int genes, string arch){//NO FUNCIONA AUN

	string File_name = "3GRN_1FB_" + arch + ".txt";
	ofstream archivo(File_name, ios::app);
	int num_pos = 0;
	int num_neg = 0;
	for(int i = 0; i < db_size; i++){
		get_nwgrn_original(nombreDB, i, genes);
		num_pos = 0;
		num_neg = 0;

		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				if(grn.nw[z][j] == 1 && j == z){
				num_pos++;
				}
			}
		}
		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				if(grn.nw[z][j] == -1 && j == z){
				num_neg++;
				}
			}
		}
		archivo << num_pos << "\t"
				<< num_neg
				<< endl;

	}
	archivo.close();
	
}

void Space::get_2fb_clusters(string nombreDB, int db_size, int genes, string arch){

	string File_name = "4GRN_2FB_" + arch + ".txt";
	ofstream archivo(File_name, ios::app);
	int num = 0;
	for(int i = 0; i < db_size; i++){
		get_nwgrn_original(nombreDB, i, genes);
		num = 0;

		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				if(grn.nw[z][j] == 1 && grn.nw[j][z] == 1 && j != z){
				num++;
				for(int k = 0; k < genes; k++){
					if(k != j && k != z){
						get_numb_conec_extra(grn.nw, genes, k);
					}
				}
				get_auto_int(grn.nw, genes, j, z);
				}
			}
		}
		archivo << num  << "\t"
		<< pos_extra_int << "\t"
		<< neg_extra_int << "\t"
		<< pos_extra_get << "\t"
		<< neg_extra_get << "\t"
		<< auto_pos_extra << "\t"
		<< auto_neg_extra << "\t"
		<< auto_pos_fb << "\t"
		<< auto_neg_fb << endl;

	}
	archivo.close();
	
}

void Space::get_3fb_clusters(string nombreDB, int db_size, int genes, string arch){

	string File_name = "3GRN_3FB_" + arch + ".txt";
	ofstream archivo(File_name, ios::app);
	int num = 0;
	for(int i = 0; i < db_size; i++){
		get_nwgrn_original(nombreDB, i, genes);
		num = 0;

		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				if(grn.nw[z][j] == 1 && j != z){
				for(int k = 0; k < genes; k++){
					if(grn.nw[k][z] == 1 && k != j && k != z){
					if(grn.nw[j][k] == 1){
						num++;
					}
					}
				}
				}
			}
		}
		if(num > 0){
			cout << i << endl;
		}
		archivo << num << endl;


	}
	archivo.close();
	
}

void Space::get_4fb_clusters(string nombreDB, int db_size, int genes){

	int num = 0;
	for(int i = 0; i < db_size; i++){
		cout << "Numero de red:" << i << endl;
		get_nwgrn_original(nombreDB, i, genes);
		grn.printnw(cout);
		num = 0;

		for(int j = 0; j < genes; j++){
			for(int z = 0; z < genes; z++){
				if(grn.nw[z][j] != 0 && j != z){	
				for(int k = 0; k < genes; k++){
					if(grn.nw[k][z] != 0 && k != j && k != z){
					for(int w = 0; w < genes; w++){
						if(grn.nw[w][k] != 0 && w != j && w != z && w != k){
						if(grn.nw[j][w] != 0){
							cout << j << z << k << w << endl;
							num++;
						}
						}
					}
					}
				}
				}
			}
		}
	}
		cout << num << endl;
	
}

void Space::get_numb_conec_extra(int **mat, int genes, int gen_extra){

	pos_extra_int = 0;
	neg_extra_int = 0;
	pos_extra_get = 0;
	neg_extra_get = 0;
	auto_pos_extra = 0;
	auto_neg_extra = 0;

	for(int i = 0; i < genes; i++){
		if(mat[i][gen_extra] < 0 && i != gen_extra){
			neg_extra_int++;
		}
		else if (mat[i][gen_extra] > 0 && i != gen_extra){
			pos_extra_int++;
		}
		
	}

	for(int i = 0; i < genes; i++){
		if(mat[gen_extra][i] < 0 && i != gen_extra){
			neg_extra_get++;
		}
		else if(mat[gen_extra][i] > 0 && i != gen_extra){
			pos_extra_get++;
		}
	}

	if(mat[gen_extra][gen_extra] == 1){
		auto_pos_extra++;
	}
	else if(mat[gen_extra][gen_extra] == -1){
		auto_neg_extra++;
	}


}

void Space::get_auto_int(int **mat, int genes, int gen1, int gen2){

	auto_pos_fb = 0;
	auto_neg_fb = 0;

	if(mat[gen1][gen1] == 1){
		auto_pos_fb++;
	}
	if(mat[gen1][gen1] == -1){
		auto_neg_fb++;
	}
	if(mat[gen2][gen2] == 1){
		auto_pos_fb++;
	}
	if(mat[gen2][gen2] == -1){
		auto_neg_fb++;
	}

}

void Space::get_list_1scc(int genes, int db_size){

	string File_name = "Num_1scc_" + bas.inttostring(genes) + ".txt";
	ofstream fa(File_name, ios::app);

	int** nw_temp = new int*[genes];
	for (int i = 0; i < genes; i++) {
    	nw_temp[i] = new int[genes];
	}
	bas.fillmatm1(nw_temp, genes, genes);

	grn.make_nw(genes, true);
	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			grn.nw[i][j] = nw_temp[i][j];
		}
	}

	for(int i = 0; i < db_size+1; i++){

		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);

		if(diff.number_of_sccs() == 1 ){
			fa << i << endl;
		}

		org2.delete_grn();
		org2.delete_diff();
		
		bas.next_in_matrix(nw_temp, genes, genes, -1, 1, 1);
		grn.make_nw(genes, true);
		for(int i = 0; i < genes; i++){
			for(int j = 0; j < genes; j++){
				grn.nw[i][j] = nw_temp[i][j];
			}
		}

	}

	for (int i = 0; i < genes; i++){
		delete [] nw_temp[i];
	}
	delete [] nw_temp;

	fa.close();

}

bool Space::get_2fb_members(int gene_target, bool* gene_list, int genes){

	bool two_nodes_cluster = false;
	two_nodes_cluster = false;
	for(int j = 0; j < genes; j++){
		for(int z = 0; z < genes; z++){
			if(grn.nw[z][j] == 1 && grn.nw[j][z] == 1 && j != z){
			two_nodes_cluster = true;
			gene_list[j] = true;
			gene_list[z] = true;
			}
		}
	}

	return two_nodes_cluster;
	
}

bool Space::get_3fb_members(int gene_target, bool* gene_list, int genes){

	bool three_nodes_cluster = false;
	three_nodes_cluster = false;

	for(int j = 0; j < genes; j++){
		for(int z = 0; z < genes; z++){
			if(grn.nw[z][j] == 1 && j != z){
			for(int k = 0; k < genes; k++){
				if(grn.nw[k][z] == 1 && k != j && k != z){
				if(grn.nw[j][k] == 1){
					three_nodes_cluster = true;
					gene_list[j] = true;
					gene_list[z] = true;
					gene_list[k] = true;
				}
				}
			}
			}
		}
	}

	return three_nodes_cluster;
	
}
void Space::conection_gen_extra(bool* gene_list, int genes){

	pos_from_extra_scc = 0;
	neg_from_extra_scc = 0;
	pos_from_scc_ext = 0;
	neg_from_scc_ext = 0;

	for(int i = 0; i < genes; i++){
	if(gene_list[i] == false){
		for(int j = 0; j < genes; j++){
		if(gene_list[j] == true){
			if(grn.nw[j][i] == 1)
				pos_from_extra_scc++;
			else if(grn.nw[j][i] == -1)
				neg_from_extra_scc++;
			else if(grn.nw[i][j] == 1)
				pos_from_scc_ext++;
			else if(grn.nw[i][j] == -1)
				neg_from_scc_ext++;
		}
		}
	}	
	}

}
bool Space::node_negative_connections_only(bool* gene_list, int genes){

	bool one_node_all_negative = false; 
	bool node_all_are_negative;

	for(int i = 0; i < genes; i++){
	if(gene_list[i] == false){
		node_all_are_negative = true;
		for(int j = 0; j < genes; j++){
		if(gene_list[j] == true){
			if(grn.nw[j][i] == 1)
				node_all_are_negative = false;
		}
		}
		if(node_all_are_negative == true){
			one_node_all_negative = true;
			return one_node_all_negative;
			break;
		}
	}
	}
}

bool Space::all_nodes_are_fb(bool* gene_list, int genes){

	int suma = 0;
	bool all_nodes_true = false;
	for(int i = 0; i < genes; i++){
		if(gene_list[i] == true)
		suma++;
	}
	if(suma == genes)
	all_nodes_true = true;

	return all_nodes_true;
}

void Space::two_three_nodes(string nombreDB, int iteration, int genes, int db_size, int part_size, string arch){
	//Quitar arch y db_size
	string File_name = "output_chunks/Analisis_" + bas.inttostring(iteration) + "_" + arch + ".txt";
	//string File_name = "archivo.txt";
	ofstream archivo(File_name, ios::app);
	ifstream file(nombreDB);

	bool* gene_list = new bool[genes];
	int columns = genes*genes;
	int i_col;

	int** nw_complete = new int*[part_size];
	for(int i = 0; i < part_size; i++){
		nw_complete[i] = new int[columns];
	}
	bas.fillmat0(nw_complete, part_size, columns);

    int row = 0;
	string line;
    while(getline(file, line) && row < part_size) {
        int col = 0;
        for (size_t i = 0; i < line.size() && col < columns; i++) {

            if(line[i] == '-'){
                nw_complete[row][col] = -1;
                col++;
                i++;
            }
            else if(line[i] == '0'){
                nw_complete[row][col] = 0;
                col++;
            }
			else if(line[i] == '1'){
                nw_complete[row][col] = 1;
                col++;
			}
        }

        row++;
    }

	grn.make_nw(genes, true);

 	for(int i = 0; i < part_size; i++){

		i_col = 0;
		for(int j = 0; j < genes; j++){
			for(int k = 0; k < genes; k++){
				grn.nw[j][k] = nw_complete[i][i_col];
				i_col++;
			}
		}
		for(int h = 0; h < genes; h++){
			gene_list[h] = false;
		}

 		archivo << get_2fb_members(i, gene_list, genes) << "\t"//i no sirve para nada
				<< get_3fb_members(i, gene_list, genes) << "\t"
				<< get_4fb_members(i, gene_list, genes) << "\t";
		rev_nodes(gene_list, genes);

		conection_gen_extra(gene_list, genes);
		archivo << pos_from_extra_scc << "\t"
				<< neg_from_extra_scc << "\t"
				<< pos_from_scc_ext << "\t"
				<< neg_from_scc_ext << "\t"
				<< node_negative_connections_only(gene_list, genes) << "\t" 
				<< all_nodes_are_fb(gene_list, genes) << endl;
	}

	archivo.close();
}

bool Space::get_4fb_members(int gene_target, bool* gene_list, int genes){

	bool four_nodes_cluster = false;
	four_nodes_cluster = false;

	for(int j = 0; j < genes; j++){
	for(int z = 0; z < genes; z++){
		if(grn.nw[z][j] == 1 && j != z){
		for(int k = 0; k < genes; k++){
			if(grn.nw[k][z] == 1 && k != j && k != z){
			for(int l = 0; l < genes; l++){
				if(grn.nw[l][k] == 1 && l != k && l != j && l != z){
				if(grn.nw[j][l] == 1){
					four_nodes_cluster = true;
					gene_list[j] = true;
					gene_list[z] = true;
					gene_list[k] = true;
					gene_list[l] = true;
				}
				}
			}
			}
		}
		}
	}
	}

	return four_nodes_cluster;
	
}

int Space::rev_nodes(bool* gene_list, int genes){

	ofstream archivo_nodo_solo("Nodo_apagado.txt", ios::app);
	int suma_entrada = 0;
	int salida_positiva = 0;
	
	for(int i = 0; i < genes; i++){
		suma_entrada = 0;
		salida_positiva = 0;
		
			for(int j = 0; j < genes; j++){
				if(grn.nw[j][i] == 1)
				salida_positiva++;
			}
			for(int j = 0; j < genes; j++){
				suma_entrada = suma_entrada + grn.nw[i][j];
			}
			archivo_nodo_solo << suma_entrada << "\t"
							  << salida_positiva << "\t";
			
		
	}
	archivo_nodo_solo << endl;

	archivo_nodo_solo.close();
}

void Space::space_distribution(string nombreDB, int genes, int db_size){

	int size = (genes*genes)+1;
	int* vector = new int[size];
	bas.fillv0(vector, size);

	int** new_nw;
    new_nw = new int*[db_size];
    for(int i = 0; i < db_size; i++){
        new_nw[i] = new int[genes*genes];
    }
    bas.fillmat0(new_nw, db_size, genes*genes);

    ifstream file(nombreDB);

	if (!file) {
    cerr << "Error: no se pudo abrir el archivo: " << nombreDB << endl;
    return;
	}
	
	int row = 0;
    string line;
    while(getline(file, line) && row < db_size){
        int col = 0;
        for (size_t i = 0; i < line.size() && col < genes*genes; i++) {

            if(line[i] == '-'){
                new_nw[row][col] = -1;
                col++;
                i++;
            }
            else if(line[i] == '0'){
                new_nw[row][col] = 0;
                col++;
            }
            else if(line[i] == '1'){
                new_nw[row][col] = 1;
                col++;
            }
	    }
        row++;
    }

	for(int i = 0; i < 20; i++){
		for(int j = 0; j < genes*genes; j++){
			cout << new_nw[i][j];
		}
		cout << endl;
	}
	cout << endl;

	int diferencias;
	for(int i = 0; i < db_size; i++){
		cout << "iteracion" << i << endl;
		//get_nwgrn_original(nombreDB, i, genes);
		for(int j = i+1; j < db_size; j++){
			//cout << "comparacion" << endl;
			//get_nwgrn_db(nombreDB, j, genes);
			//grn_check.printnw(cout);
			diferencias = 0;

			for(int z = 0; z < genes*genes; z++){
				if(new_nw[i][z] != new_nw[j][z]){
					diferencias++;
				}
			}

			vector[diferencias]++;

		}

	}

	for(int i = 0; i < size; i++){
		cout << vector[i] << endl;
	}

}

void Space::remove_cicatt(string nombreDB, string new_nombreDB, int genes, int db_size){

	ofstream fa(new_nombreDB, ios::app);

	for(int i = 0; i < db_size; i++){
		get_nwgrn_original(nombreDB, i, genes);
		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		if(org2.cicatt == 0){
			diff.printnw(cout);	
			for(int j = 0; j < genes; j++){
				for(int x = 0; x < genes; x++){
					fa << grn.nw[j][x];
				}
			}
			fa << endl;
		}
		org2.delete_grn();
		org2.delete_diff();	
	}

	fa.close();

}

void Space::seek_cicatt_allspace(int db_size, int genes){
	
	//	Recorre todo el espacio y busca si hay atractores ciclicos, si los hay agraga al contador
	//	funciona por que el primero todos -1 nunca tiene ciclos.

	grn.make_nw(genes, true);
	for(int i = 0; i < genes; i++){
		for(int j = 0; j < genes; j++){
			grn.fill_nw_value(-1, i, j);
		}
	}
	grn.printnw(cout);

	numero_targets = 1; 

	for(int i = 1; i < db_size; i++){
		cout << "Numero de iteracion: " << i << endl;
		bas.next_in_matrix(grn.nw, genes, genes, -1, 1, 1);
		org2.make_org(grn, genes);
		diff.copy(org2.diff_org);
		if(org2.cicatt == 0  && diff.number_of_sccs() == 1 ){
			numero_targets++;
		}
		org2.delete_grn();
		org2.delete_diff();

	}
	cout << numero_targets << endl;

}

void Space::make_graphviz(string nombreDB, int genes, int num_red){

    //grn.make_nw(genes, true);
	for(int i = 0; i < 1497; i++){
		get_nwgrn_original(nombreDB, i, genes);
		
		grn.print_dot_edgecolor("Figuras/regulacion" + bas.inttostring(i));
		bas.get_dot_fig("Figuras/regulacion" + bas.inttostring(i));
	}

}

void Space::genes_comp(string nombreDB, string nombreDB2, int db_size, int genes){

	int numero_positivas_ab = 0; 
	int numero_negativas_ab = 0;
	int numero_positivas_ba = 0; 
	int numero_negativas_ba = 0;
	int pos_mismo_grupo = 0;

	string File_name = "Num_pos_neg.txt";
	ofstream fa(File_name, ios::app);

	for(int x = 0; x < db_size; x++){
		get_group(nombreDB, x);
		get_nwgrn_original(nombreDB2, x, genes);

		numero_positivas_ab = 0; 
		numero_negativas_ab = 0;
		numero_positivas_ba = 0; 
		numero_negativas_ba = 0;
		pos_mismo_grupo = 0;


		for(int i = 0; i < grupo_1.size(); i++){
			for(int j = 0; j < grupo_2.size(); j++){
				if(grn.nw[grupo_2[j]][grupo_1[i]] == 1)
				numero_positivas_ab++;
				else if(grn.nw[grupo_2[j]][grupo_1[i]] == -1)
				numero_negativas_ab++;
			}
		}

		for(int i = 0; i < grupo_2.size(); i++){
			for(int j = 0; j < grupo_1.size(); j++){
				if(grn.nw[grupo_1[j]][grupo_2[i]] == 1)
				numero_positivas_ba++;
				else if(grn.nw[grupo_1[j]][grupo_2[i]] == -1)
				numero_negativas_ba++;
			}
		}

		for(int i = 0; i < grupo_1.size(); i++){
			for(int j = 0; j < grupo_1.size(); j++){
				if(i != j){
					if((grn.nw[grupo_1[j]][grupo_1[i]] == 1) && (grn.nw[grupo_1[i]][grupo_1[j]] == 1)){
						pos_mismo_grupo++;
					}
				}
			}
		}

		for(int i = 0; i < grupo_2.size(); i++){
			for(int j = 0; j < grupo_2.size(); j++){
				if(i != j){
					if((grn.nw[grupo_2[j]][grupo_2[i]] == 1) && (grn.nw[grupo_2[i]][grupo_2[j]] == 1)){
						pos_mismo_grupo++;
					}
				}
			}
		}

		if(pos_mismo_grupo == 0){
			for(int i = 0; i < grupo_1_original.size(); i++){
				for(int j = 0; j < grupo_1_original.size(); j++){
					if(i != j){
						if((grn.nw[grupo_1_original[j]][grupo_1_original[i]] == 1) && (grn.nw[grupo_1_original[i]][grupo_1_original[j]] == 1)){
							pos_mismo_grupo++;
						}
					}
				}
			}
			for(int i = 0; i < grupo_2_original.size(); i++){
				for(int j = 0; j < grupo_2_original.size(); j++){
					if(i != j){
						if((grn.nw[grupo_2_original[j]][grupo_2_original[i]] == 1) && (grn.nw[grupo_2_original[i]][grupo_2_original[j]] == 1)){
							pos_mismo_grupo++;
						}
					}
				}
			}
		}

		cout << "Matriz:" << x << endl;
		for (int i = 0; i < grupo_1.size(); i++) {
    		cout << grupo_1[i] << " ";
		}
		cout << endl;
		for (int i = 0; i < grupo_2.size(); i++) {
    		cout << grupo_2[i] << " ";
		}
		cout << endl;

		fa << numero_positivas_ab << "\t"
		<< numero_negativas_ab << "\t"
		<< numero_positivas_ba << "\t"
		<< numero_negativas_ba << "\t"
		<< pos_mismo_grupo << endl;

		cout << numero_positivas_ab << "\t"
		<< numero_negativas_ab << "\t"
		<< numero_positivas_ba << "\t"
		<< numero_negativas_ba << "\t"
		<< pos_mismo_grupo << endl;
		
	}
	fa.close();
}

void Space::get_group(string nombreDB, int target){

    ifstream archivo(nombreDB);
    string linea;

    for (int i = 0; i <= target; i++) {
        getline(archivo, linea);
    }
    grupo_1.clear();
    grupo_2.clear();

    size_t tab = linea.find('\t');
    string columna_1 = linea.substr(0, tab);
    string columna_2 = linea.substr(tab + 1);

    istringstream valores_1(columna_1);
    int valor, indice = 0;

    while (valores_1 >> valor) {
        if (valor == 1)
            grupo_1.push_back(indice);
        indice++;
    }

    istringstream valores_2(columna_2);
    indice = 0;

    while (valores_2 >> valor) {
        if (valor == 1)
            grupo_2.push_back(indice);
        indice++;
    }

	grupo_1_original.clear();
	grupo_2_original.clear();
	grupo_1_original = grupo_1;
	grupo_2_original = grupo_2;

	vector<int> grupo_1_filtrado;
	vector<int> grupo_2_filtrado;

	for (int valor : grupo_1) {
		if (find(grupo_2.begin(), grupo_2.end(), valor) == grupo_2.end()) {
			grupo_1_filtrado.push_back(valor);
		}
	}
	for (int valor : grupo_2) {
		if (find(grupo_1.begin(), grupo_1.end(), valor) == grupo_1.end()) {
			grupo_2_filtrado.push_back(valor);
		}
	}

	grupo_1 = grupo_1_filtrado;
	grupo_2 = grupo_2_filtrado;

}
