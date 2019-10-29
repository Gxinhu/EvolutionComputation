package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF14 extends DF {
	private double a, b, G, Wt, Nt, Ht, kt, pt;
	double[] y = new double[variables];

	public DF14(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 3, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF14";
		lowerLimit_[0] = 0;
		upperLimit_[0] = 1;
		lowerLimit_[1] = 0;
		upperLimit_[1] = 1;
		for (int var = 2; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 2; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - G, 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.sin(0.5 * Math.PI * t);
		double[] f = new double[objectives];
		y[0] = 0.5 + G * (x.getValue(0) - 0.5);
		f[0] = gx(x) * (1 - y[0] + 0.05 * Math.sin(6 * Math.PI * y[0]));
		f[1] = gx(x) * (1 - x.getValue(1) + 0.05 * Math.sin(6 * Math.PI * x.getValue(1)))
				* (y[0] + 0.05 * Math.sin(6 * Math.PI * y[0]));
		f[2] = gx(x) * (x.getValue(1) + 0.05 * Math.sin(6 * Math.PI * x.getValue(1)))
				* (y[0] + 0.05 * Math.sin(6 * Math.PI * y[0]));
		solution.setObjective(0, f[0]);
		solution.setObjective(1, f[1]);
		solution.setObjective(2, f[2]);
	}

	@Override
	public void dynamicChange(int iteration) {
		super.dynamicChange(iteration);
	}
}