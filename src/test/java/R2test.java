import jmetal.qualityIndicator.R2;

public class R2test {
	public static void main(String[] args) {
		double[][] archive = {{0, 1}, {0.2, 0.8}, {1, 0}, {1.8, 0.5}, {0.4, 1.2}, {1.5, 1.5}};
		double[][] lambdas = {{1, 0.001}, {0.5, 0.5}, {0.001, 1}};
		double[] R2incdicators = new double[archive.length];
		R2 r2 = new R2(2, lambdas);
		double r = r2.R2(archive, archive);
		for (int i = 0; i < archive.length; i++) {
			R2incdicators[i] = r - r2.R2Withouth(archive, archive, i);
		}
		System.out.println(R2incdicators);
	}

}
