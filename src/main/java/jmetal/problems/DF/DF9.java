package jmetal.problems.DF;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.archive.DominateArchive;
import jmetal.util.wrapper.XReal;

public class DF9 extends DF {
	private double a, b, G, Wt, Nt;
	double[] y = new double[variables];

	public DF9(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DF9";
		lowerLimit_[0] = 0;
		upperLimit_[0] = 1;
		for (int var = 1; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -1;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 1.0;
		for (int i = 1; i < variables; i++) {
			sums += Math.pow(x.getValue(i) - Math.cos(4 * t + x.getValue(0) + x.getValue(i - 1)), 2);
		}
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		Nt = 1 + Math.floor(10 * Math.abs(Math.sin(0.5 * t * Math.PI)));
		double[] f = new double[objectives];
		f[0] = gx(x) * (x.getValue(0) + Math.max(0, (1 / (2 * Nt) + 0.1) * Math.sin(2 * Nt * Math.PI * x.getValue(0))));
		f[1] = gx(x) * (1 - x.getValue(0) + Math.max(0, (1 / (2 * Nt) + 0.1) * Math.sin(2 * Nt * Math.PI * x.getValue(0))));
		solution.setObjective(0, f[0]);
		solution.setObjective(1, f[1]);
	}

	@Override
	public void dynamicChange(int iteration) {
		super.dynamicChange(iteration);
	}

	@Override
	public double[][] getPF() {
		SolutionSet solutions = new SolutionSet(numOfPF2d);
		double[] x = new double[numOfPF2d];
		for (int i = 0; i < numOfPF2d; i++) {
			x[i] = (double) i / (numOfPF2d - 1);
			solutions.add(new Solution(objectives));
		}
		double[][] f = new double[numOfPF2d][objectives];
		for (int i = 0; i < numOfPF2d; i++) {
			for (int j = 0; j < objectives; j++) {
				if (j == 0) {
					f[i][j] = x[i] + Math.max(0, (0.1 + 0.5 / Nt) * Math.sin(2 * Nt * Math.PI * x[i]));
					solutions.get(i).setObjective(j, f[i][j]);
				} else if (j == 1) {
					f[i][j] = 1 - x[i] + Math.max(0, (0.1 + 0.5 / Nt) * Math.sin(2 * Nt * Math.PI * x[i]));
					solutions.get(i).setObjective(j, f[i][j]);
				}
			}
		}
		DominateArchive archive = new DominateArchive(numOfPF2d, objectives);
		for (int i = 0; i < solutions.size(); i++) {
			archive.add(solutions.get(i));
		}
		return archive.writeObjectivesToMatrix();
	}
}