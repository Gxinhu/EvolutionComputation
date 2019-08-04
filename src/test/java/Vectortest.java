import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;

public class Vectortest {
	public static void main(String[] args) {
		double[][] lambda_;
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(4, 1, 6);
		lambda_ = vg.getVectors();
		System.out.println(lambda_[0][0]);
	}
}
