package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.archive.DominateArchive;
import jmetal.util.wrapper.XReal;

public class DF12 extends DF {
	private double a, b, G, Wt, Nt, Ht, kt;
	double[] y = new double[variables];

	public DF12(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 3, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF12";
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
			sums += Math.pow(x.getValue(i) - t * x.getValue(0), 2);
		}
		sums += Math.abs(Math.sin(Math.floor(kt * (2 * x.getValue(0) - 1)) * Math.PI / 2) * Math.sin(Math.floor(kt * (2 * x.getValue(1) - 1)) * Math.PI / 2));
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		kt = 10 * Math.sin(Math.PI * t);
		double[] f = new double[objectives];
		f[1] = gx(x) * Math.pow(Math.sin(0.5 * Math.PI * x.getValue(1)) * Math.cos(0.5 * Math.PI * x.getValue(0)), 1);
		f[0] = gx(x) * Math.pow(Math.cos(0.5 * Math.PI * x.getValue(1)) * Math.cos(0.5 * Math.PI * x.getValue(0)), 1);
		f[2] = gx(x) * Math.pow(Math.sin(0.5 * Math.PI * x.getValue(0)), 1);
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
		SolutionSet solutions = new SolutionSet(numOfPF * numOfPF);
		double[][] x = new double[2][numOfPF];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < numOfPF; j++) {
				x[i][j] = (double) j / (numOfPF - 1);
			}
		}
		double[][] temp = new double[numOfPF][numOfPF];
		for (int i = 0; i < numOfPF; i++) {
			for (int j = 0; j < numOfPF; j++) {
				temp[i][j] = Math.abs(Math.sin(Math.floor(kt * (2 * x[0][i] - 1)) * Math.PI / 2) * Math.sin(Math.floor(kt * (2 * x[1][j] - 1)) * Math.PI / 2)) + 1;
				solutions.add(new Solution(objectives));
			}
		}
		double[][] f = new double[numOfPF * numOfPF][objectives];
		for (int i = 0; i < numOfPF; i++) {
			for (int k = 0; k < numOfPF; k++) {
				for (int j = 0; j < objectives; j++) {
					if (j == 0) {
						f[i * numOfPF + k][j] = temp[i][k] * Math.cos(0.5 * Math.PI * x[1][k]) * Math.cos(0.5 * Math.PI * x[0][i]);
						solutions.get(i * numOfPF + k).setObjective(j, f[i * numOfPF + k][j]);
					} else if (j == 1) {
						f[i * numOfPF + k][j] = temp[i][k] * Math.sin(0.5 * Math.PI * x[1][k]) * Math.cos(0.5 * Math.PI * x[0][i]);
						solutions.get(i * numOfPF + k).setObjective(j, f[i * numOfPF + k][j]);
					} else if (j == 2) {
						f[i * numOfPF + k][j] = temp[i][k] * Math.sin(0.5 * Math.PI * x[0][i]);
						solutions.get(i * numOfPF + k).setObjective(j, f[i * numOfPF + k][j]);
					}
				}
			}
		}
		DominateArchive archive = new DominateArchive(numOfPF * numOfPF, objectives);
		for (int i = 0; i < solutions.size(); i++) {
			archive.add(solutions.get(i));
		}
		return archive.writeObjectivesToMatrix();
	}
}