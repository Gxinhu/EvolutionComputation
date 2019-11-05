package jmetal.problems.dynamicProblem;

import jmetal.core.Problem;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.PseudoRandom;

public abstract class dynamicProblem extends Problem {

	protected int variables;
	protected int objectives;
	protected int severityOfChanges;
	protected int numberOfChanges;
	protected int t0;
	protected double t;
	protected int numOfPF2d = 1500;
	protected int numOfPF = 50;
	protected int r = 0;

	public dynamicProblem(String solutionType, int varables, int objctives, int severrityOfchanges, int numberOfChanges, int t0) {
		this.variables = varables;
		this.objectives = objctives;
		this.severityOfChanges = severrityOfchanges;
		this.numberOfChanges = numberOfChanges;
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
		int temp = Math.max(iteration + numberOfChanges - (t0), 0);
		double lastT = t;
		t = ((double) 1 / severityOfChanges) * Math.floor((double) temp / numberOfChanges);
		if (t != lastT) {
			r = PseudoRandom.randInt(0, variables - 1);
		}
	}

	@Override
	public double[][] getPF() {
		double[][] PF = {{0}, {0}};
		return PF;
	}
}
