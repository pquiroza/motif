#include<cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
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


int generateRandomnumber(int limite){

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0,limite);

  //cout << "Generando " << uni(rng) << '\n';
  return uni(rng);
   //static const double fraction = 1.0 / (RAND_MAX + 1.0);
   //return 0 + static_cast<int>((limite - 0 + 1) * (std::rand() * fraction));

}


Matrix generaMatrizInicial(vector<string> datos, int largo){
  //cout << datos.size() << '\n';
  vector<Kmer> kmers;
  for (int i = 0;i<datos.size();i++){


    int t = generateRandomnumber(datos.at(i).length()-largo);


    string pal = datos.at(i).substr(t,largo);
    Kmer k = Kmer(pal,t,largo);
    kmers.push_back(k);

  }

Matrix m = Matrix(kmers,0,"");
return m;

}

void showMatriz(Matrix m){
  for (int i =0;i<m.kmers.size();i++){
    cout << m.kmers.at(i).adn<< " " << m.kmers.at(i).posicioninicial << '\n';
  }
}



/*
template<typename K, typename V>
void print_map(std::unordered_map<K,V> const &m)
{
    for (auto const& pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}
*/

float getFitness(Matrix m){

  unordered_map<char,int> umap;

int suma = 0;

for (int i =0;i<m.kmers.at(0).adn.length();i++){

  umap['A'] = 0;
  umap['C'] = 0;
  umap['E'] = 0;
  umap['D'] = 0;
  umap['G'] = 0;
  umap['F'] = 0;
  umap['I'] = 0;
  umap['H'] = 0;
  umap['K'] = 0;
  umap['M'] = 0;
  umap['L'] = 0;
  umap['N'] = 0;
  umap['Q'] = 0;
  umap['P'] = 0;
  umap['S'] = 0;
  umap['R'] = 0;
  umap['T'] = 0;
  umap['W'] = 0;
  umap['V'] = 0;
  umap['Y'] = 0;
  umap['X'] = 0;

for (int l =0; l<m.kmers.size();l++){
  umap[m.kmers.at(l).adn[i]] = umap[m.kmers.at(l).adn[i]] + 1;
}
int maximo = 0;
int ceros = 0;
for (auto& it: umap) {
    // Do stuff
    if (it.second>maximo){
      maximo = it.second;
    }
    if (it.second==0){
      ceros++;
    }
  //  cout << it.first << " " << it.second << '\n';
}
suma = suma + maximo +  ceros;






}
float total = (float) suma / (m.kmers.at(0).adn.length() * (m.kmers.size()+(m.kmers.size()-1)));
return total;
}

void genetico(vector<string> datos, int poblacion, int ciclos, int lkmer,string salida){
  vector<Matrix> pop;
  for (int i =0;i<poblacion;i++) {
  Matrix m = generaMatrizInicial(datos,lkmer);
  pop.push_back(m);
}


for (int i =0;i<pop.size();i++){
  cout << "Fitness " << getFitness(pop.at(i)) << '\n';
}

}




int main(int nNumberofArgs, char* pszArgs[])
{



vector<string> datos =cargaDatos("PS00010.fa");
//cout << "Datos : " << datos.size() << '\n';

genetico(datos,100000,100,6,"salida");
/*
for (int i=0; i<datos.size();i++){
  cout << datos.at(i) << '\n';

}
*/


}
