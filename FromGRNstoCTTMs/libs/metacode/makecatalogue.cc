#include <fstream>
#include <iostream>
#include <time.h>
#include <math.h>
#include<stdio.h>
#include "/Users/ce/libs/alea.h"
#include "/Users/ce/libs/basics.h"
#include "/Users/ce/libs/graphi.h"
#include "/Users/ce/libs/graphb.h"
#include "/Users/ce/libs/graphc.h"


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

void construye(int tam, bool dirono);
void set_pars(int tam, bool dirono, int& nupre, int& previos);
int get_props_nw(GraphB red, int tam, bool dirono, int* mochi, int previos);

int main(int argc, char **argv)
{//previos corregir
	
	//directed


	construye(2, true);//
	construye(2, false);
	construye(3, false);
//	construye(3, true);
//
	construye(4, false);


  return 0;
}

void set_pars(int tam, bool dirono, int& nupre, int& previos)
{
	if (dirono) {
		if (tam==2) {
			nupre = 7;
			previos = 2;
		}
		if (tam==3) {
			nupre = 86;
			previos = 9;
		}
		if (tam==4) {
			nupre = 1780; //?
			previos = 95;
		}
		
	}
	else {
		if (tam==2) {
			nupre = 3;
			previos = 2;
		}
		if (tam==3) {
			nupre = 10;
			previos = 5;
		}
		if (tam==4) {
			nupre = 50; 
			previos = 15;
		}
		
	}

	
}

int get_props_nw(GraphB red, int tam, bool dirono, int* mochi, int previos)
{ //previos es necesario??
	int i;
	
	int ed = red.number_of_edges();
	mochi[0] = red.count_non_selfinteracting_nodes();
	mochi[1] = red.count_selfinteracting_nodes();
	if (dirono) {
		if (tam >= 3) {
			red.count_subgraphs_2();
			for (i=2; i<9; i++) {
				mochi[i] = red.how_many_subgraphs(2, (i-2));
			}			
		}
		if (tam == 4) {
			red.count_subgraphs_3();
			for (i=9; i<95; i++) {
				mochi[i] = red.how_many_subgraphs(3, (i-9));
			}
		}
			
	}
	else {
		if (tam >= 3) {
			red.count_subgraphs_2();
			for (i=2; i<5; i++) {
				mochi[i] = red.how_many_subgraphs(2, (i-2));
			}
		}
		if (tam == 4) {
			red.count_subgraphs_3();
			for (i=5; i<15; i++) {
				mochi[i] = red.how_many_subgraphs(3, (i-5));
			}
		}
		
		
	}

	

	return ed;
}

void construye(int tam, bool dirono)
{
  Alea jacta(1);
  jacta.check_rng(cout);
	Basics bas(jacta);
  int i, j, nupre,van,com,previos, ede, k;
  GraphB red(jacta);
	int *cont;
	GraphB *limo;
	bool ya;
	int *mochi, **limochi, *vedge;
	ofstream fs, fz, fc;
	string cu,er,da,oo,so,ga;
	bool **matriz;
	set_pars(tam, dirono, nupre, previos);//
	vedge = new int[nupre];
	bas.fillv0(vedge, nupre);
	cont = new int[nupre];
	bas.fillv0(cont, nupre);
	mochi = new int[previos];
	bas.fillv0(mochi, previos);
	limochi = new int*[nupre];
	if (limochi==NULL) 
		bas.escape("no alcanza");
	for (i=0; i<nupre; i++) {
		limochi[i] = new int[previos];
		if (limochi[i]==NULL)
			bas.escape("no alcanza");
	}
	bas.fillmatm1(limochi, nupre, previos);
	limo = new GraphB[nupre];
	if (limo==NULL)
		bas.escape("no alcanzo.");
	for (i=0; i<nupre; i++) 
		limo[i].start_rng(jacta);
	van =0;
	matriz = new bool*[tam];
	if (matriz==NULL) 
		bas.escape("no llego");
	for (i=0; i<tam; i++) {
		matriz[i] = new bool[tam];
		if (matriz[i]==NULL)
			bas.escape("no llego");
	}
	bas.fillmat0(matriz, tam, tam);
	
	
	do {
		if (dirono)
			bas.next_in_matrix(matriz, tam, tam);
		else {
			do {
				bas.next_in_matrix(matriz, tam, tam);
			} while (!bas.is_symmetric(matriz, tam, tam));
		}
		red.make_nw(tam, dirono);
		if (dirono) {
			for (i=0; i<tam; i++) 
				for (j=0; j<tam; j++) 
					red.force_interaction(j, i, matriz[i][j]);
		}
		else {
			for (i=0; i<tam; i++) 
				for (j=0; j<=i; j++) 
					red.force_interaction_undir(j, i, matriz[i][j]);
		}
		com = red.get_components();
		if (com==1) {
			ya = false;
			ede = get_props_nw(red, tam, dirono, mochi, previos);
			for (i=0; i<van; i++) {
				if (ede == vedge[i]) {
					if (bas.eqvec(mochi, previos, limochi[i], previos)) {
						if (limo[i].isomorphic_pair(red)) {
							cont[i]++;		
							ya = true;
							break;
						}
					}
				}
			}
			if (!ya) {
				limo[van].copy(red);
				vedge[van]=ede;
				cont[van]=1;
				for (j=0; j<previos; j++) 
					limochi[van][j] = mochi[j];
				van++; 	
			}
		}
		red.clear();
		cout << van << endl;
	} while (!bas.last_matrix(matriz, tam, tam));
		
		
	if (dirono)
		cu = "Directed/";
	else 
		cu = "Undirected/";
	er = bas.inttostring(tam);
	da = cu+er+".tex";
	oo = cu+er+".cc";
	
	
	if (nupre != van) {
		cout << "[Error]: AH JIJO she dijo\n";
		exit(1);
	}
	
	
	fz.open(da.c_str()); //tex
//	fs.open(oo.c_str()); // prog
	
	
	
	fz << "\\documentclass[11pt]{article}\n\\usepackage{graphicx}\n\\usepackage[utf8]{inputenc}\n\\usepackage{url}\n\\usepackage{times}\n\\usepackage{placeins}\n\\usepackage{amssymb,amsmath}\n\\usepackage[usenames]{color}\n\\usepackage{multirow}\n\\usepackage{array}\n\\usepackage{colortbl}\n";
	
	fz << "\\textheight 9.3in\n\\textwidth 6.3in\n";
	fz << "\\begin{document}\n";
	fz << "\\begin{table}[!h]\n\\centering\n\\begin{tabular}[!h]{|c| c| c| c|}\n";
	fz <<"\\hline\n";
	fz << "Id&Subgraph&edges&count\\\\\\hline\n";
	
//	string nombre, cunupre, nompremot, numnompremot, socont;
//
//	nombre = "LM_"+er;
//	
//	cunupre = "nupre_"+er;
//	nompremot = "nompremot_"+er;
//	numnompremot ="numnompremot_"+er;
//	socont = "socont_"+er;
//	
//
//	fs << cunupre << "=" << nupre << ";\n";
//	fs << nombre << "= new GraphB[" << cunupre << "];\n";
//	fs << "if (" << nombre << "==NULL)\n\tbasic.escape(\"a graph vector\");\n";
//	fs << "for (i=0; i<" << cunupre << "; i++) {\n";
//	fs << "\t" << nombre << "[i].start_rng(est);\n";
//	fs << "\t" << nombre << "[i].make_nw(" << tam << ", ";
//	if (dirono)
//		fs << "true";
//	else 
//		fs << "false";
//	fs << ");\n}\n";
//	
//
//	fs << nompremot << "= new int*[" << cunupre << "];\n";
//	fs << "if (" << nompremot << "==NULL)\n\tbasic.escape(\"a graph vector\");\n";
//	fs <<  "for (i=0; i<" << cunupre << "; i++) {\n";
//	fs << "\t" << nompremot << "[i]= new int[" << numnompremot << "];\n";
//	fs << "\tif (" << nompremot << "[i]==NULL)\n\tbasic.escape(\"a graph vector\");\n}\n";
//
//	fs << socont << "= new int[" << cunupre << "];\n";
//	fs << "if (" << socont << "==NULL)\n\tbasic.escape(\"a graph vector\");\n";
//
//	fs << "basic.fillv0(" << socont << ", " << cunupre << ");\n";
//	fs << "basic.fillmat0(" << nompremot << ", " << cunupre << ", " << numnompremot << ");\n";
//	
	for (i=0; i<nupre; i++) {
//		fs << socont << "[" << i << "] = " << cont[i] << ";\n";
//		for (j=0; j< previos; j++) 
//			fs << nompremot << "[" << i << "][" << j << "] = " << limochi[i][j] << ";\n";
//		if (dirono) {
//			for (j=0; j<tam; j++) 
//				for (k=0; k<tam; k++) 
//					if (limo[i].weight(k, j)) 
//						fs << nombre << "[" << i << "].force_interaction(" << k << ", " << j << ", true);\n";
//		}
//		else {
//			for (j=0; j<tam; j++) 
//				for (k=0; k<=j; k++) 
//					if (limo[i].weight(k, j)) 
//						fs << nombre << "[" << i << "].force_interaction_undir(" << k << ", " << j << ", true);\n";
//		}
		so = bas.inttostring(i);
		ga = cu+er+"_"+so;
		limo[i].print_dot_wo_labels_circ(ga, 1.0);
		bas.get_dot_fig(ga);
		fz << i << "&\\mbox{\\includegraphics [width=0.08\\textwidth]{" << cu << tam << "_" << i << ".pdf}}&" << vedge[i] << "&" << cont[i] << "\\\\\\hline\n";
		if (((i+1)%13)==0) {
			fz << "\\end{tabular}\n\\end{table}\n\\newpage\n\n";
			
			fz << "\\begin{table}[!h]\n\\centering\n\\begin{tabular}[!h]{|c| c| c| c|}\n";
			fz <<"\\hline\n";
			fz << "Id&Subgraph&edges&count\\\\\\hline\n";//estaba por aqui
		}
	}
	fz << "\\end{tabular}\n\\end{table}\n\\newpage\n\n\\end{document}\n";
	fz.close();
//	fs.close(); //faltaba
	string orden;
	orden = "/usr/local/texlive/2011/bin/x86_64-darwin/pdflatex -output-directory ";
	orden = orden+cu+" "+er+".tex";
	system(orden.c_str());
	
	orden = "rm "+cu+"*.log "+cu+"*.aux";
	
	system(orden.c_str());
	
}
