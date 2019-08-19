package motif.classes;
import java.io.BufferedReader;
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

public class GeneticMotif {






  static ArrayList<String> datos = new ArrayList<String>();

    static String linea ="";
    public static int largo = 6;
    public static int poblacion = 100000;
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



  public static Matriz generaInicial(){
    int indices[] = new int[datos.size()];
    for (int i=0;i<datos.size();i++){
      Random r = new Random();
      int indice = r.nextInt(datos.get(i).length()-largo);
      indices[i] = indice;

    }
    Matriz m = new Matriz(indices,0);
    return m;
  }


  public static double getFitness(Matriz m){
    String[] palabras = new String[datos.size()];
    int parcial = 0;


    for (int i =0;i<m.indices.length;i++){
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
  return propor;
  }


public static ArrayList<Matriz> roulette(ArrayList<Matriz> pop){
  ArrayList<Matriz> ganadores = new ArrayList<Matriz>();

  for (int i=0;i<pop.size();i++){
    Random r = new Random();
    double compara = r.nextDouble();
    if (pop.get(i).fitness>compara){
      ganadores.add(pop.get(i));
    }

  }
  System.out.println("ganadores "+ganadores.size());
  return ganadores;

}

public static ArrayList<Matriz> cruzamiento(Matriz m1,Matriz m2){
  int[] temporal = new int[datos.size()];
  int[] temporal2 = new int[datos.size()];

  Random r = new Random();
  int indice = r.nextInt(datos.size());

for (int i=0;i<m1.indices.length;i++){
  System.out.println(m1.indices[i]+" "+m2.indices[i]);
}
System.out.println("-----------------------");

  for(int i=0;i<indice;i++){
    temporal[i] = m1.indices[i];
    temporal2[i] = m2.indices[i];
  }
  for (int l=indice;l<datos.size();l++){
    temporal[l]=m2.indices[l];
    temporal2[l] = m1.indices[l];
  }
for (int i=0;i<temporal.length;i++){
  System.out.println(temporal[i]+" "+temporal2[i]);
}
return null;

}

  public static void main(String args[]){







    ArrayList<Matriz> pop = new ArrayList<Matriz>();
    System.out.println("Java");
    cargaDatos("../PS00821.fa");


    for (int i=0;i<poblacion;i++){
    Matriz m = generaInicial();
    m.fitness = getFitness(m);

    pop.add(m);
  }

pop = roulette(pop);
for (int i=0;i<pop.size();i++){
  Random r = new Random();
  int indice = r.nextInt(pop.size());
  int indice2 = r.nextInt(pop.size());

  cruzamiento(pop.get(indice),pop.get(indice2));


}

System.out.println(pop.size());



  }
}
