#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphb.h"
#include "graphi.h"
#include "graphi0.h"
#include "Organismo2.h"
#include "Spaces.h"

using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv){

    Organismo_2 org2;
    GraphI0 grn;
    Basics bas;
    Space spa;


    int** nw_complete = new int*[2];//Tamaño de la matriz este archivo solo tenia 2
    for(int i = 0; i < 2; i++){
        nw_complete[i] = new int[144];//Numero de entradas de la matriz para 12 = 144
    }
    bas.fillmat0(nw_complete, 2, 144);

    for(int iteration = 1; iteration < 2; iteration++){
        //string File_name = "output_chunks/chunk_" + bas.inttostring(iteration) + ".txt";
        string File_name = "12random_nw_100mil_con oreja.txt";
        ifstream file(File_name);
        int i_col;

        int row = 0;
        string line;
        while(getline(file, line) && row < 2){
            int col = 0;
            for (size_t i = 0; i < line.size() && col < 144; i++) {

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
        grn.make_nw(12, true);
        for(int i = 0; i < 2; i++){
            i_col = 0;
            for(int j = 0; j < 12; j++){
                for(int k = 0; k < 12; k++){
                    grn.nw[j][k] = nw_complete[i][i_col];
                    i_col++;
                }
            }
        grn.printnw(cout);
        org2.make_multi_org(grn, 12, iteration);
        } 
    }

    for(int i = 0; i < 2; i++){
        delete nw_complete[i];
    }
    delete [] nw_complete;

    return 0;
}