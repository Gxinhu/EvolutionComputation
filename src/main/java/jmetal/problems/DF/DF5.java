package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF5 extends DF {
	private double a, b, G, Wt;

	public DF5(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF5";
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
			sums += Math.pow(x.getValue(i) - G, 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.sin(0.5 * t * Math.PI);
		Wt = Math.floor(10 * G);
		double[] f = new double[objectives];
		f[0] = gx(x) * (x.getValue(0) + 0.02 * Math.sin(Wt * Math.PI * x.getValue(0)));
		f[1] = gx(x) * (1 - x.getValue(0) + 0.02 * Math.sin(Wt * Math.PI * x.getValue(0)));
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
					f[i][j] = x[i] + 0.02 * Math.sin(Wt * Math.PI * x[i]);
				} else if (j == 1) {
					f[i][j] = 1 - x[i] + 0.02 * Math.sin(Wt * Math.PI * x[i]);
				}
			}
		}
		return f;
	}
}