package jmetal.problems.singleobjective.cec2005;

public class getcec2005bias {
	public double bias;
	private int number;

	public getcec2005bias(int number) {
		this.number = number;
		this.bias = getcec2005biass();
	}

	public double getcec2005biass() {
		double[] bias = new double[]{
				-4.5000000e+002, -4.5000000e+002, -4.5000000e+002, -4.5000000e+002, -3.1000000e+002,
				3.9000000e+002, -1.8000000e+002, -1.4000000e+002, -3.3000000e+002, -3.3000000e+002,
				9.0000000e+001, -4.6000000e+002, -1.3000000e+002, -3.0000000e+002, 1.2000000e+002, 1.2000000e+002,
				1.2000000e+002, 1.0000000e+001, 1.0000000e+001, 1.0000000e+001, 3.6000000e+002, 3.6000000e+002, 3.6000000e+002,
				2.6000000e+002, 2.6000000e+002
		};
		return bias[number];
	}
}
