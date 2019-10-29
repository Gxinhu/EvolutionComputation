package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DF11 extends DF {
	private double a, b, G, Wt, Nt, Ht;
	double[] y = new double[variables];

	public DF11(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 3, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF11";
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 1.0 + G;
		for (int i = 2; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - 0.5 * G * x.getValue(0), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.abs(Math.sin(0.5 * Math.PI * t));
		for (int i = 0; i < 2; i++) {
			y[i] = G * Math.PI / 6 + (Math.PI / 2 - G * Math.PI / 3) * x.getValue(i);
		}
		double[] f = new double[objectives];
		f[0] = gx(x) * Math.sin(y[0]);
		f[1] = gx(x) * Math.sin(y[1]) * Math.cos(y[0]);
		f[2] = gx(x) * Math.cos(y[0]) * Math.cos(y[1]);
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
		double[][] y = new double[2][numOfPF];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < numOfPF; j++) {
				y[i][j] = Math.PI * G / 6 + (Math.PI / 2 - Math.PI * G / 3) * x[i][j];
			}
		}
		double[][] f = new double[numOfPF * numOfPF][objectives];
		for (int i = 0; i < numOfPF; i++) {
			for (int k = 0; k < numOfPF; k++) {
				for (int j = 0; j < objectives; j++) {
					if (j == 0) {
						f[i * numOfPF + k][j] = Math.sin(y[0][i]);
					} else if (j == 1) {
						f[i * numOfPF + k][j] = Math.sin(y[1][k]) * Math.cos(y[0][i]);
					} else if (j == 2) {
						f[i * numOfPF + k][j] = Math.cos(y[1][k]) * Math.cos(y[0][i]);
					}
				}
			}
		}
		return f;
	}
}