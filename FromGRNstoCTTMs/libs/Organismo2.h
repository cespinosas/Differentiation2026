#ifndef ORGANISMO2_H
#define ORGANISMO2_H

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

class Organismo_2
{
public:

    Organismo_2();
    Organismo_2(Alea& seed);

    // ---- Hacer red diferenciacion ---- //
    void star_rng(Alea& seed);

    void make_org(GraphI0& grn, int genes);
    void make_multi_org(GraphI0& grn, int genes, int arch);
	void initial_vector(int genes);
	void contador();
	void num_hileras();
	void new_find_all_att(int genes);
	void new_fill_diff_cell_nw(int genes);

    void delete_grn();
    void delete_diff();

    void get_all_att(string Archivo, int genes);
    void make_org_zero(GraphI0& grn2, int genes, int itera, int arch);
    void zero_att_position(string Archivo, int genes, int val, int num_nw);

    GraphI0 grn_org;
	GraphC diff_org;
    GraphC diff_org_check;
    int fil;
    int cicatt;

private:

    Alea ale;
    Basics bas;

    int eracciones;
    int elev;
	//int cicatt;
	int hileras;

	bool sies;
	bool yaesta;

    int* ve;
	int* cont;
	int* phil;
	int** provatt;
	int*** allatt;

};

#endif