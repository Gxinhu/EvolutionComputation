package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF3 extends DF {
	public DF3(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF3";
		lowerLimit_[0] = 0;
		lowerLimit_[1] = 1;
		for (int var = 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1;
			upperLimit_[var] = 2;
		}
	}

	public double G() {
		return Math.sin(0.5 * Math.PI * t);
	}

	public double H() {
		return 1.5 + G();
	}

	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 1; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - G() - Math.pow(x.getValue(0), H()), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		double[] f = new double[objectives];
		f[0] = x.getValue(0);
		f[1] = gx(x) * (1 - Math.pow(x.getValue(0) / gx(x), H()));
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
					f[i][j] = 1 - Math.pow(x[i], H());
				}
			}
		}
		return f;
	}
}