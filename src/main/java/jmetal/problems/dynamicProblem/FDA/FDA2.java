package jmetal.problems.dynamicProblem.FDA;

import jmetal.core.Solution;
import jmetal.problems.dynamicProblem.dynamicProblem;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class FDA2 extends dynamicProblem {
	private double G, Ht;

	public FDA2(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "FDA2";
		lowerLimit_[0] = 0;
		upperLimit_[0] = 1;
		for (int var = 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 1; i < variables - 7; i++) {
			sums += Math.pow(x.getValue(i), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		Ht = 2 * Math.sin(0.5 * Math.PI * (t - 1));
		double[] f = new double[objectives];
		f[0] = x.getValue(0);
		double sumss = Ht;
		for (int i = variables - 7; i < variables; i++) {
			sumss += Math.pow(x.getValue(i) - Ht / 4, 2);
		}
		f[1] = gx(x) * (1 - Math.pow(x.getValue(0) / gx(x), Math.pow(2, sumss)));
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
					f[i][j] = 1 - Math.pow(f[i][j - 1], Math.pow(2, Ht));
				}
			}
		}
		return f;
	}
}