package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF6 extends DF {
	private double a, b, G, Wt;
	double[] y = new double[variables];

	public DF6(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF6";
		lowerLimit_[0] = 0;
		upperLimit_[0] = 1;
		for (int var = 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1;
			upperLimit_[var] = 1;
		}
	}

	public double H() {
		return 1.5 + a;
	}

	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 1; i < variables; i++) {
			sums += Math.abs(G) * Math.pow(y[i], 2) - 10 * Math.cos(2 * Math.PI * y[i]) + 10;
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.sin(0.5 * t * Math.PI);
		for (int i = 0; i < variables; i++) {
			y[i] = x.getValue(i) - G;
		}
		a = 0.2 + 2.8 * Math.abs(G);
		double[] f = new double[objectives];
		f[0] = gx(x) * Math.pow(x.getValue(0) + 0.1 * Math.sin(3 * Math.PI * x.getValue(0)), a);
		f[1] = gx(x) * Math.pow(1 - x.getValue(0) + 0.1 * Math.sin(3 * Math.PI * x.getValue(0)), a);
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
					f[i][j] = Math.pow(x[i] + 0.1 * Math.sin(3 * Math.PI * x[i]), a);
				} else if (j == 1) {
					f[i][j] = Math.pow(1 - x[i] + 0.1 * Math.sin(3 * Math.PI * x[i]), a);

				}
			}
		}
		return f;
	}
}