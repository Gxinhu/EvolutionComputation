package jmetal.problems.dynamicProblem.FDA;

import jmetal.core.Solution;
import jmetal.problems.dynamicProblem.dynamicProblem;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class FDA4 extends dynamicProblem {
	private double G, Ht, Ft;

	public FDA4(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 3, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "FDA4";

		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 0;
		for (int i = 2; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - G, 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		G = Math.abs(Math.sin(0.5 * Math.PI * t));
		double[] f = new double[objectives];
		f[0] = (1 + gx(x)) * Math.cos(0.5 * Math.PI * x.getValue(0)) * Math.cos(0.5 * Math.PI * x.getValue(1));
		f[1] = (1 + gx(x)) * Math.cos(0.5 * Math.PI * x.getValue(0)) * Math.sin(0.5 * Math.PI * x.getValue(1));
		f[2] = (1 + gx(x)) * Math.sin(Math.PI * 0.5 * x.getValue(0));
		solution.setObjective(0, f[0]);
		solution.setObjective(1, f[1]);
		solution.setObjective(2, f[2]);
	}

	@Override
	public void dynamicChange(int iteration) {
		super.dynamicChange(iteration);
	}

	//	@Override
	public double[][] getPsF() {
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
					f[i][j] = (1 + G) * (1 - Math.pow(f[i][j - 1] / (1 + G), 0.5));

				}
			}
		}
		return f;
	}

	@Override
	public double[][] getPF() {
		double[][] x = new double[2][numOfPF];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < numOfPF; j++) {
				x[i][j] = (double) j * (Math.PI / 2) / (numOfPF - 1);
			}
		}

		double[][] f = new double[numOfPF * numOfPF][objectives];
		for (int i = 0; i < numOfPF; i++) {
			for (int k = 0; k < numOfPF; k++) {
				for (int j = 0; j < objectives; j++) {
					if (j == 0) {
						f[i * numOfPF + k][j] = Math.cos(x[0][i]) * Math.cos(x[1][k]);
					} else if (j == 1) {
						f[i * numOfPF + k][j] = Math.cos(x[0][i]) * Math.sin(x[1][k]);
					} else if (j == 2) {
						f[i * numOfPF + k][j] = Math.sin(x[0][i]);
					}
				}
			}
		}
		return f;
	}
}