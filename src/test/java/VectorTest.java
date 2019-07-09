import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class VectorTest {
	public static void main(String[] args) {
		double[] temp = {1, 2, 3};
		RealVector vector = new ArrayRealVector(temp);
		vector = vector.mapAdd(0.1);
		System.out.println(vector.toArray());
	}
}
