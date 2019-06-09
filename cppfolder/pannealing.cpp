#include<cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


class Kmer {
public:
  string adn;
  int posicioninicial;
  int largo;
};

class Matrix
{ public:
  vector<Kmer> kmers;
  double fitnness;
  string code;

};

vector<string> cargaDatos(string archivo){
  string line;
  vector<string> datos;
  string secuencia ="";
  ifstream archivofa (archivo);
  if (archivofa.is_open()){
    while(getline(archivofa,line)){

      if (line.at(0)=='>'){
        if(secuencia==""){

        }
        else {
          cout << "Insetando datos" << '\n';
          datos.push_back(secuencia);
            secuencia = "";
        }

      }
      else{
        //cout << line << '\n';
        secuencia = secuencia + line;

      }

    }

    archivofa.close();
  }
  else {
    cout << "no se puede abrir archivo";
  }
  cout << "Datos : " << datos.size();
  for (int i=0; i<datos.size();i++){
    cout << datos.at(i) << '\n';
  }

  return datos;



}


int main(int nNumberofArgs, char* pszArgs[])
{
cargaDatos("PS00010.fa");



}
