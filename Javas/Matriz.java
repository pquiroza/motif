package motif.classes;
import java.lang.Cloneable;
public class Matriz implements Cloneable{
  public int[] indices;
  public double fitness;


public Matriz(int[] indices,double fitness){
  this.indices = indices;
  this.fitness = fitness;

}

public  Object clone() throws CloneNotSupportedException {
return super.clone();
}

}
