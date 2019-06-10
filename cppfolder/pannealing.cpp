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

  Kmer(string a, int pi, int l){
    adn = a;
    posicioninicial = pi;
    largo = l;
  }
};

class Matrix
{ public:
  vector<Kmer> kmers;
  float fitness;
  string code;

  Matrix(vector<Kmer> k,float fts, string cod){
    kmers = k;
    fitness = fts;
    code = cod;
  }


  void setKmer(int posicion,Kmer kmer){
    kmers.at(posicion) = kmer;
  }
  void setFitness(float fitness2) {
    fitness = fitness2;
  }

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


  return datos;



}


Matrix generaMatrizInicial(vector<string> datos, int largo){

}


int main(int nNumberofArgs, char* pszArgs[])
{
vector<string> datos =cargaDatos("PS00010.fa");
cout << "Datos : " << datos.size();
/*
for (int i=0; i<datos.size();i++){
  cout << datos.at(i) << '\n';

}
*/


}
