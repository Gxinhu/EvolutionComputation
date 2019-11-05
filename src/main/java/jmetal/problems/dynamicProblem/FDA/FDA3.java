package jmetal.problems.dynamicProblem.FDA;

import jmetal.core.Solution;
import jmetal.problems.dynamicProblem.dynamicProblem;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class FDA3 extends dynamicProblem {
	private double G, Ht, Ft;

	public FDA3(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "FDA3";
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
		double sums = 1.0 + G;
		for (int i = 2; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - G, 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		Ft = Math.pow(10, 2 * Math.sin(0.5 * Math.PI * t));
		G = Math.abs(Math.sin(0.5 * Math.PI * t));
		double[] f = new double[objectives];
		f[0] = (Math.pow(x.getValue(0), Ft) + Math.pow(x.getValue(1), Ft)) / 2;
		f[1] = gx(x) * (1 - Math.pow(f[0] / gx(x), 0.5));
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
					f[i][j] = (1 + G) * (1 - Math.pow(f[i][j - 1] / (1 + G), 0.5));
				}
			}
		}
		return f;
	}

	//	@Override
	public double[][] getsPF() {
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
						f[i * numOfPF + k][j] = (Math.pow(x[0][i], Ft) + Math.pow(x[1][k], Ft)) / 2;
					} else if (j == 1) {
						f[i * numOfPF + k][j] = (1 + G) * (1 - Math.pow(f[i * numOfPF + k][j - 1] / (1 + G), 0.5));
					}
				}
			}
		}
		return f;
	}
}