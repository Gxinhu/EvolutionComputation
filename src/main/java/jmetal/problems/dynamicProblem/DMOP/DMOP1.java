package jmetal.problems.dynamicProblem.DMOP;

import jmetal.core.Solution;
import jmetal.problems.dynamicProblem.dynamicProblem;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class DMOP1 extends dynamicProblem {
	private double G, Ht;

	//这里的上下界我看论文之间有区别 这里的测速问题是按照
	//A Steady-State and Generational Evolutionary Algorithm for Dynamic Multiobjective Optimization
	public DMOP1(String solutionType, int varables, int severrityOfchanges, int numberOfChanges, int t0) {
		super(solutionType, varables, 2, severrityOfchanges, numberOfChanges, t0);
		problemName_ = "DMOP1";
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = 0;
			upperLimit_[var] = 1;
		}
	}


	public double gx(XReal x) throws JMException {
		double sums = 0;
		for (int i = 1; i < variables; i++) {
			sums += Math.pow(x.getValue(i), 2);
		}
		sums += 1 + (9 / variables) * sums;
		return sums;
	}

	@Override
	public void evaluate(Solution solution) throws JMException {
		XReal x = new XReal(solution);
		Ht = 0.75 * Math.sin(0.5 * Math.PI * t) + 1.25;
		double[] f = new double[objectives];
		f[0] = x.getValue(0);
		f[1] = gx(x) * (1 - Math.pow(x.getValue(0) / gx(x), Ht));
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
					f[i][j] = 1 - Math.pow(x[i], Ht);
				}
			}
		}
		return f;
	}
}