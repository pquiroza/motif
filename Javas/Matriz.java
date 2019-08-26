package motif.classes;
import java.lang.Cloneable;
public class Matriz implements Cloneable{
  public int[] indices;
  public double fitness;
  public int state;


public Matriz(int[] indices,double fitness, int state){
  this.indices = indices;
  this.fitness = fitness;
  this.state = state;

}

public  Object clone() throws CloneNotSupportedException {
return super.clone();
}

}
