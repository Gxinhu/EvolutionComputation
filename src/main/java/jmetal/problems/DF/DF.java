package jmetal.problems.DF;

import jmetal.core.Problem;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;

public abstract class DF extends Problem {

	protected int variables;
	protected int objectives;
	protected int numberOfChange;
	protected int serverityOfChange;
	protected int t0;
	protected double t;
	protected int numOfPF2d = 1500;
	protected int numOfPF = 50;

	public DF(String solutionType, int varables, int objctives, int serverityOfChange, int numberOfChange, int t0) {
		this.variables = varables;
		this.objectives = objctives;
		this.numberOfChange = numberOfChange;
		this.serverityOfChange = serverityOfChange;
		this.t0 = t0;
		setNumberOfVariables(varables);
		numberOfObjectives_ = this.objectives;
		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		if (solutionType.compareTo("BinaryReal") == 0) {
			solutionType_ = new BinaryRealSolutionType(this);
		} else if (solutionType.compareTo("Real") == 0) {
			solutionType_ = new RealSolutionType(this);
		} else {
			System.out.println("Error: solution type " + solutionType + " invalid");
			System.exit(-1);
		}

	}

	@Override
	public void dynamicChange(int iteration) {
		int temp = Math.max(iteration + numberOfChange - (t0), 0);
		t = ((double) 1 / serverityOfChange) * Math.floor((double) temp / numberOfChange);
	}

	@Override
	public double[][] getPF() {
		double[][] PF = {{0}, {0}};
		return PF;
	}
}