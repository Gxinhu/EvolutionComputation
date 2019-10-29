package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF10 extends DF {
	private double a, b, G, Wt, Nt, Ht;
	double[] y = new double[variables];

	public DF10(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 3, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF10";
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
			sums += Math.pow(x.getValue(i) - Math.sin(2 * Math.PI * (x.getValue(0) + x.getValue(1))) / (1 + G), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		Ht = 2.25 + 2 * Math.cos(0.5 * Math.PI * t);
		G = Math.sin(0.5 * Math.PI * t);
		double[] f = new double[objectives];
		f[0] = gx(x) * Math.pow(Math.sin(0.5 * Math.PI * x.getValue(0)), Ht);
		f[1] = gx(x) * Math.pow(Math.sin(0.5 * Math.PI * x.getValue(1)) * Math.cos(0.5 * Math.PI * x.getValue(0)), Ht);
		f[2] = gx(x) * Math.pow(Math.cos(0.5 * Math.PI * x.getValue(1)) * Math.cos(0.5 * Math.PI * x.getValue(0)), Ht);
		solution.setObjective(0, f[0]);
		solution.setObjective(1, f[1]);
		solution.setObjective(2, f[2]);
	}

	@Override
	public void dynamicChange(int iteration) {
		super.dynamicChange(iteration);
	}

	@Override
	public double[][] getPF() {
		double[][] x = new double[2][numOfPF];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < numOfPF; j++) {
				x[i][j] = (double) j / (numOfPF - 1);
			}
		}
		double[][] f = new double[numOfPF * numOfPF][objectives];
		for (int i = 0; i < numOfPF; i++) {
			for (int k = 0; k < numOfPF; k++) {
				for (int j = 0; j < objectives; j++) {
					if (j == 0) {
						f[i * numOfPF + k][j] = Math.pow(Math.sin(x[0][i] * 0.5 * Math.PI), Ht);
					} else if (j == 1) {
						f[i * numOfPF + k][j] = Math.pow(Math.sin(x[1][k] * 0.5 * Math.PI), Ht) * Math.pow(Math.cos(x[0][i] * 0.5 * Math.PI), Ht);
					} else if (j == 2) {
						f[i * numOfPF + k][j] = Math.pow(Math.cos(x[1][k] * 0.5 * Math.PI), Ht) * Math.pow(Math.cos(x[0][i] * 0.5 * Math.PI), Ht);

					}
				}
			}
		}
		return f;
	}
}