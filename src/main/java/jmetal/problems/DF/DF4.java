package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF4 extends DF {
	private double a, b;

	public DF4(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF4";
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -2;
			upperLimit_[var] = 2;
		}
	}

	public double H() {
		return 1.5 + a;
	}

	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 1; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - (a * Math.pow(x.getValue(0), 2) / i), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		a = Math.sin(0.5 * Math.PI * t);
		b = 1 + Math.abs(Math.cos(0.5 * Math.PI * t));
		double[] f = new double[objectives];
		f[0] = gx(x) * Math.pow(Math.abs(x.getValue(0) - a), H());
		f[1] = gx(x) * Math.pow(Math.abs(x.getValue(0) - a - b), H());
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
			x[i] = a + (double) (b) * i / (numOfPF2d - 1);
		}
		double[][] f = new double[numOfPF2d][objectives];
		for (int i = 0; i < numOfPF2d; i++) {
			for (int j = 0; j < objectives; j++) {
				if (j == 0) {
					f[i][j] = Math.pow(Math.abs(x[i] - a), H());
				} else if (j == 1) {
					f[i][j] = Math.pow(b - Math.pow(f[i][j - 1], 1 / H()), H());

				}
			}
		}
		return f;
	}
}