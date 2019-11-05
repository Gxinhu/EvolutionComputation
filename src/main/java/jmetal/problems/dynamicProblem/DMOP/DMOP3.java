package jmetal.problems.dynamicProblem.DMOP;

import jmetal.core.Solution;
import jmetal.problems.dynamicProblem.dynamicProblem;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DMOP3 extends dynamicProblem {
	private double G, Ht;

	public DMOP3(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DMOP3";
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 1;
		for (int i = 0; i < variables; i++) {
			if (i != r) {
				sums += Math.pow(x.getValue(i) - G, 2);
			}
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.abs(Math.sin(0.5 * Math.PI * t));
		double[] f = new double[objectives];
		f[0] = x.getValue(r);
		f[1] = gx(x) * (1 - Math.pow(x.getValue(r) / gx(x), 0.5));
		solution.setObjective(0, f[0]);
		solution.setObjective(1, f[1]);
	}

	@Override
	public void dynamicChange(int iteration) {
		super.dynamicChange(iteration);
	}

	@Override
	public double[][] getPF() {
		double[] x = new double[numOfPF2d];
		for (int i = 0; i < numOfPF2d; i++) {
			x[i] = (double) i / (numOfPF2d - 1);
		}
		double[][] f = new double[numOfPF2d][objectives];
		for (int i = 0; i < numOfPF2d; i++) {
			for (int j = 0; j < objectives; j++) {
				if (j == 0) {
					f[i][j] = x[i];
				} else if (j == 1) {
					f[i][j] = 1 - Math.pow(x[i], 0.5);
				}
			}
		}
		return f;
	}

}