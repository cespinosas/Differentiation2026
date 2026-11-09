#ifndef SPACE_H
#define SPACE_H

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
#include "graphc.h"
#include "Organismo2.h"
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm> 
#include <utility>


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

class Space
{
public:

    Space();
    Space(Alea& seed);
    void star_rng(Alea& seed);

    int return_grn(int i, int j);
    int return_matriz_nodo();
    int return_nw_positive_interactions();
    void get_n_nw(string Archivo, int filas);
    void get_positive_interactions(string Archivo, int i);
    void get_nwgrn_db(string Archivo, int i, int genes);//lo guarda en grn_check
    void get_nwgrn_original(string Archivo, int i, int genes);//lo guarda en grn
    void save_grn(string File_name, int genes);

    void make_all_nw(const string& nombresal, const string& nombre_pos_intera, const string& nombre_nodo_intera, const string& nombrentr, int n, int start, int end);
    void check_one_component();
    void check_interaction_rows(int fil_col);
    void check_nodes_interaction(int fil_col);
    void nw_combination(int fil_col);
    void func_combination(int fil_col);
    void check_nwdb_interactions(const string& nombreDB, int fil_col, int id);

    //----- Comparar BD vs estructura pre establecida ----
    
    void review_tar_nw(string nombreDB, int db_size, int genes);
    void get_distance_db(string nombreDB, int db_size, int genes, int att);
    void comp_all_space(string nombreDB, int genes);
    void ofunc_comb(int fil);
    void nw_combination_dif(int** nw, int fil_col);

    //----- Comparar dos bases de datos en txt ------
    void comp_two_db(string nombreDB1, string nombreDB2, string nombresal, int size_db1, int size_db2, int fil_col);
    bool del_1row(string nombreDB, int genes);

    //----- Reduce el tamaño de una BD partiendola en partes iguales ------
    void red_db(string nombreDB, int part_size, int genes, int itera);

    //----- Busca el tamaño de los pozos en el espacio ------
    void make_random_wells();
    void space_wells_txt(string nombreDB, string nombreDB2, int db_size, int genes);
    void space_wells(string nombreDB, string nombreDB2, int db_size, int genes);

    void random_walks(string nombreDB, int genes, int itera);
    void random_walks_2ears(string nombreDB, int genes, int itera);
    void random_nw(string nombreDB, int genes, int itera);

    void change_all_well(int **mat, int well_number_seek, int new_well_number, int db_size);
    void get_nwgrn_memory(int **mat, int i, int genes);
    void get_nwgrncheck_memory(int **mat, int i, int genes);
    void review_mutation(int **mat, int db_size, int genes);
    void fill_false(string nombreDB, string nombreDB2);
    void change_state(string nombreDB, int target_line);
    void random_mutation(int genes);
    void check_mutation(string nombreDB, int db_size, int genes);
    void get_next_idx(string nombreDB, int db_size);
    void count_one_step(string nombreDB, string nombreDB2, int genes, int db_size1, int db_size2);
    void number_of_diff(int genes);
    void count_number_differences(string nombreDB, string nombreDB2, int genes, int db_size1, int db_size2);


    //-----Busca feedback loops de distitntos tamaños en las redes-------
    void get_1fb_clusters(string nombreDB, int db_size, int genes, string arch);
    void get_2fb_clusters(string nombreDB, int db_size, int genes, string arch);
    void get_3fb_clusters(string nombreDB, int db_size, int genes, string arch);
    void get_4fb_clusters(string nombreDB, int db_size, int genes);//AUN SE NECESITA REVISAR

    void get_numb_conec_extra(int **mat, int genes, int gen_extra);
    void get_auto_int(int **mat, int genes, int gen1, int gen2);

    //-----Buscar SCC en las redes------
    void get_list_1scc(int genes, int db_size);

    bool get_2fb_members(int gene_target, bool* gene_list, int genes);
    bool get_3fb_members(int gene_target, bool* gene_list, int genes);
    bool get_4fb_members(int gene_target, bool* gene_list, int genes);
    void conection_gen_extra(bool* gene_list, int genes);
    bool node_negative_connections_only(bool* gene_list, int genes);
    bool all_nodes_are_fb(bool* gene_list, int genes);
    void two_three_nodes(string nombreDB, int iteration, int genes, int db_size, int part_size, string arch);
    int rev_nodes(bool* gene_list, int genes);

    void space_distribution(string nombreDB, int genes, int db_size);
    void remove_cicatt(string nombreDB, string new_nombreDB, int genes, int db_size);
    void seek_cicatt_allspace(int db_size, int genes);
    void make_graphviz(string nombreDB, int genes, int num_red);

    void genes_comp(string nombreDB, string nombreDB2, int db_size, int genes);
    void get_group(string nombreDB, int target);

    GraphI0 grn;
    GraphC diff;


private:

    Alea ale;
    Basics bas;
	GraphI0 grn_check;
	GraphC conden;
    GraphC diff_check;
    GraphC diff_ear;
    Organismo_2 org2;

    int genes;
    int elev;
    int num_interactions_newnw;
    int num_interactions_bdnw;
    int num_positive_interactions_newnw;
    int num_positive_interactions_db;
    int numero_targets;
    int saved_value;
    int row1;
    int row2;
    int idx_nw;
    //int label_target;
    int idx_target;
    int well_size;
    int flag;
    int number_diff;
    int pos_extra_int;
    int neg_extra_int;
    int pos_extra_get;
    int neg_extra_get;
    int auto_pos_extra;
    int auto_neg_extra;
    int auto_pos_fb;
    int auto_neg_fb;

    int pos_from_scc_ext;
    int pos_from_extra_scc;
    int neg_from_scc_ext;
    int neg_from_extra_scc;
    
    bool yaparams;
    bool save_nw;
    bool nw_target;
    bool one_component;
    bool check_decrease_rows;
    bool no_empty_rows;
    bool one_step;

    int** nw_combinations;

    vector<int> grupo_1, grupo_2;
    vector<int> grupo_1_original, grupo_2_original;

};

#endif