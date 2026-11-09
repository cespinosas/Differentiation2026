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
#include "Spaces.h"
#include "Organismo2.h"

using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv){

  int seed = atoi(argv[1]);

  Alea cion(seed);
  Basics bas(cion);
  GraphC grac;
  GraphI0 graio;
  Organismo_2 org;
  Space spa(cion);

  string nombre = "Random_walk_6_100milpasos_sin_oreja_sin_repetidos.txt";
  string nombre2 = "1.txt";
  //spa.random_walks(nombre, 6, 100000);//Caminata aleatoria a partir de una red y comprueba si son el target
  //spa.random_nw(nombre, 6, 100000);//Hace x numero de redes aleatorias y comprueba si son el target
  //spa.random_walks_2ears(nombre, 5, 100000);
  //spa.space_distribution(nombre, 6, 39866);//Compara todas con todas las redes para ver que tan lejos estan
  //spa.remove_cicatt(nombre, nombre2, 3, 5827);//toma una BD y le quita todos los ciclos, luego devuelve la base limpia sin ciclos
  spa.seek_cicatt_allspace(81, 2);//Busca atractores ciclicos en todo el espacio y te da el numero exacto
  //19683
  //43046721


  return 0;
  
}