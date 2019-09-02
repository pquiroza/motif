package motif.classes;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;
import java.util.ArrayList;
import java.util.Random;
import motif.classes.Matriz;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.Collections;
import java.lang.Math;
public class GeneticMotif implements Cloneable{






static ArrayList<String> datos = new ArrayList<String>();
static ArrayList<Matriz> pop = new ArrayList<Matriz>();
static ArrayList<Matriz> pop2 = new ArrayList<Matriz>();


    static String linea ="";
    public static int largo = 4;
    public static int poblacion = 5000;
  public static void cargaDatos(String filename){
    try {
      Stream<String> lines = Files.lines(Paths.get(filename));
      System.out.println("Leyendo el archivo");

      lines.forEach(line->{

        if(line.startsWith(">")) {
          if (linea!=""){
            datos.add(linea);
          }
          linea="";
        }
        else{
          linea = linea+line;
        }

      });

      lines.close();

    }
    catch(IOException io){
      io.printStackTrace();
    }
  }

public static void escribeFitness(double fitness,String archivo) throws IOException {
  BufferedWriter writer = new BufferedWriter(new FileWriter(archivo,true));
  writer.write(fitness+"");
  writer.write("\n");
  writer.close();

}

public static void escribeIndividuo(Matriz m, String archivo) throws IOException{
  BufferedWriter writer = new BufferedWriter(new FileWriter(archivo));
  for(int i=0;i<m.indices.length;i++){
    writer.write(datos.get(i).substring(m.indices[i],m.indices[i]+largo)+"  "+m.indices[i]);
    writer.write("\n");
  }
  writer.write(m.fitness+"");
writer.close();
}


  public static Matriz generaInicial(){
    int indices[] = new int[datos.size()];
    for (int i=0;i<datos.size();i++){
      Random r = new Random();
      int indice = r.nextInt(datos.get(i).length()-largo);
      indices[i] = indice;

    }
    Matriz m = new Matriz(indices,0,0);
    return m;
  }


  public static void getFitness(Matriz m){
    if(m==null){
      System.out.println("ES NULL");
      System.exit(0);
    }
    String[] palabras = new String[datos.size()];
    int parcial = 0;

    //System.out.println(m.indices.length+" "+m.state);
    for (int i =0;i<datos.size();i++){

      String pal = datos.get(i).substring(m.indices[i],m.indices[i]+largo);
      palabras[i] = pal;

    }

  for (int i=0;i<largo;i++){
    HashMap<String,Integer> diccionario = new HashMap<String,Integer>();
    diccionario.put("A",0);
    diccionario.put("C",0);
    diccionario.put("D",0);
    diccionario.put("E",0);
    diccionario.put("F",0);
    diccionario.put("G",0);
    diccionario.put("H",0);
    diccionario.put("I",0);
    diccionario.put("K",0);
    diccionario.put("L",0);
    diccionario.put("M",0);
    diccionario.put("N",0);
    diccionario.put("P",0);
    diccionario.put("Q",0);
    diccionario.put("R",0);
    diccionario.put("S",0);
    diccionario.put("T",0);
    diccionario.put("V",0);
    diccionario.put("W",0);
    diccionario.put("X",0);
    diccionario.put("Y",0);
    for (int l=0;l<palabras.length;l++){
      int suma = diccionario.get(String.valueOf(palabras[l].charAt(i)))+1;
      diccionario.replace(String.valueOf(palabras[l].charAt(i)),suma);

    }


    int maximo = Collections.max(diccionario.values());
    int ceros = Collections.frequency(diccionario.values(),0);
    if (maximo == datos.size()){
      parcial = parcial + largo;
    }

    parcial = parcial + maximo +  ceros;



  }
  double propor = (double) parcial / ((datos.size()+20+largo)*largo);
  m.fitness=propor;
  //System.out.println(m.fitness);
  //return propor;
  }


public static void roulette(Matriz m){


    Random r = new Random();
    double compara = r.nextDouble();
    if (m.fitness>compara){
      m.state = 1;
    //  pop2.add(m);
    }


  //System.out.println("ganadores "+pop.size());
  //return ganadores;

}

public static void backto(Matriz m){
    m.state=0;
}

public static void cruzamiento(Matriz m1){
  int[] temporal = new int[datos.size()];
  int[] temporal2 = new int[datos.size()];

  Random r = new Random();
  int indi = r.nextInt(pop.size());
  int indice = r.nextInt(datos.size());


  for(int i=0;i<indice;i++){
    temporal[i] = m1.indices[i];
    temporal2[i] = pop.get(indi).indices[i];
  }
  for (int l=indice;l<datos.size();l++){
    temporal[l]=pop.get(indi).indices[l];
    temporal2[l] = m1.indices[l];
  }

Matriz n1 = new Matriz(temporal,0,2);
Matriz n2 = new Matriz(temporal2,0,2);

pop2.add(n1);
pop2.add(n2);

}


public static Matriz MutacionAnnealing(Matriz m){
  Random r = new Random();
  int cuantos = r.nextInt(datos.size());


  int[] nuevo = new int[datos.size()];
  for (int i=0;i<datos.size();i++){
    nuevo[i] = m.indices[i];
  }
    for (int c=0;c<cuantos;c++){
      int gen = r.nextInt(datos.size());
  int cromosoma = r.nextInt(datos.get(gen).length()-largo);
  nuevo[gen] = cromosoma;
}


  Matriz n = new Matriz(nuevo,0,0);

  return n;
}


public static void Mutacion(){

  Random r = new Random();

  //int amutar = r.nextInt((Integer) pop.size()*(0.2));
int amutar = 1000;

  for (int i=0;i<amutar;i++){

  int indi = r.nextInt(pop.size());
  int gen = r.nextInt(datos.size());

  int cromosoma = r.nextInt(datos.get(gen).length()-largo);
  pop.get(i).indices[gen] = cromosoma;
}
}


public static void GeneticRun(){
  double mejorf = 0;
  int[] psa =new int[datos.size()];
  Matriz mejori = new Matriz(psa,0,0);
   double sumatotal = 0;

  System.out.println("Java");
  String archivodatos = "PS00023.fa";
  cargaDatos("../"+archivodatos);


  for (int i=0;i<poblacion;i++){
  Matriz m = generaInicial();
  if (m==null){
    System.out.println("error inicial");
    System.exit(0);
  }


  pop.add(m);
}
pop.stream().forEach(m -> getFitness(m));
for (int c=0;c<10000000;c++){
mejorf=0;


pop.parallelStream().forEach(m-> roulette(m));
pop.removeIf(m-> m.state==0);




/*
pop.clear();
for(int l=0;l<pop2.size();l++){
Matriz mc = (Matriz) pop2.get(l).clone();
pop.add(mc);
}
pop2.clear();
*/

while(pop2.size()<poblacion){
pop.parallelStream().forEach(m -> cruzamiento(m));
}






int contanull = 0;
for(int e=0;e<pop2.size();e++){
if(pop2.get(e)!=null){
  contanull++;
  pop.add(pop2.get(e));
}
}






pop.removeIf(m-> m.state==1);
pop2.clear();



Random r = new Random();
double muta = r.nextDouble();

if(muta>0.8){
Mutacion();
}


pop.parallelStream().forEach(m -> getFitness(m));




for (int i=0;i<pop.size();i++){

sumatotal = sumatotal + pop.get(i).fitness;
if (pop.get(i).fitness > mejorf){
  mejorf = pop.get(i).fitness;


}

if(pop.get(i).fitness>mejori.fitness){
  mejori = pop.get(i);
  System.out.println("CICLO "+c+" PROMEDIO "+(double) sumatotal/pop.size()+" MEJOR "+mejorf +" POP "+pop.size());

  try {
  escribeIndividuo(mejori,"Results"+archivodatos+" "+largo+".txt");
  if(mejori.fitness==1.0){
    System.out.println("Perfect Motif Find");
    System.exit(0);
  }
}
catch(IOException e){
System.out.println("error");
}

}

}


pop2.clear();
try {
escribeFitness(mejorf,"Fitness"+archivodatos+" "+largo+".txt");

}
catch(IOException e){
System.out.println("error");
}
//System.out.println(pop.size()+" "+sumatotal);
pop.parallelStream().forEach(m -> backto(m) );
sumatotal=0;


}
}


public static void AnnealingRun(){
  String archivodatos = "PS00023.fa";
  cargaDatos("../"+archivodatos);
  double tinicial = 10;
  Matriz m = generaInicial();
  getFitness(m);
  while (tinicial>0.001){
    try {
          escribeFitness(m.fitness,"Annealing"+archivodatos+" "+largo+".txt");
        }
        catch(IOException e){

        }
    for (int c=0;c<1000;c++){


      Matriz n = MutacionAnnealing(m);
      getFitness(n);
        System.out.println("MEJOR FITNESS "+m.fitness+" "+tinicial+" "+n.fitness);
      double resta = m.fitness - n.fitness;
      //System.out.println(resta);
      if (resta<0){
        //System.out.println("ELige mejor");
        m = n;
      }
      else {
        double valor = Math.exp(-(resta/tinicial));
        Random r = new Random();
        double compara = r.nextDouble();
        if (compara<valor){
            //System.out.println("ELige peor");
          m = n;
        }

      }

    }
    tinicial = tinicial*.999;
  }
}


  public static void main(String args[]) throws CloneNotSupportedException{

AnnealingRun();



}
}
