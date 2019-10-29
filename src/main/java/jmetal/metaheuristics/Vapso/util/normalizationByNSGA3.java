
package jmetal.metaheuristics.Vapso.util;

import Jama.Matrix;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;

/**
 * This class defined a balanceable fitness estimation method(as
 * used in NMPSO).
 */
public class normalizationByNSGA3 {
	private int numberOfObjectives;
	public double[] idealPoint, nadirPoint, intercepts;
	private double[][] extremePoints;

	public normalizationByNSGA3(int numberOfObjectives) {
		this.numberOfObjectives = numberOfObjectives;
		idealPoint = new double[numberOfObjectives];
		nadirPoint = new double[numberOfObjectives];
	}

	public static void main(String[] args) {
		normalizationByNSGA3 normalization = new normalizationByNSGA3(3);
		SolutionSet front = new SolutionSet(3);
		Solution newSoution = new Solution(3);
		newSoution.setObjective(0, 0.0);
		newSoution.setObjective(1, 0.02759);
		newSoution.setObjective(2, 9.97);
		front.add(newSoution);
		Solution newSoution1 = new Solution(3);
		newSoution1.setObjective(0, 0.0);
		newSoution1.setObjective(1, 0.0);
		newSoution1.setObjective(2, 10.0);
		front.add(newSoution1);
		Solution newSoution2 = new Solution(3);
		newSoution2.setObjective(0, 0.0);
		newSoution2.setObjective(1, 0.2659);
		newSoution2.setObjective(2, 9.96);
		front.add(newSoution2);
		normalization.normalization(front);
		front.writeTransObjectivesToMatrix();
	}

	public void normalization(SolutionSet front) {
		computeIdealPoint(front);
		computeMaxPoint(front);
		computeExtremePoints(front);
		computeIntercepts();
		normalizePopulation(front);
	}

	private void computeIdealPoint(SolutionSet solutionSet) {

		for (int j = 0; j < numberOfObjectives; j++) {
			idealPoint[j] = 1.0e+30;
			for (int i = 0; i < solutionSet.size(); i++) {
				if (solutionSet.get(i).getObjective(j) < idealPoint[j]) {
					idealPoint[j] = solutionSet.get(i).getObjective(j);
				}
			}

		}

	}

	private void computeMaxPoint(SolutionSet solutionSet) {
		for (int j = 0; j < numberOfObjectives; j++) {
			nadirPoint[j] = -1.0e+30;

			for (int i = 0; i < solutionSet.size(); i++) {
				if (solutionSet.get(i).getObjective(j) > nadirPoint[j]) {
					nadirPoint[j] = solutionSet.get(i).getObjective(j);
				}
			}
		}
	}

	private void computeExtremePoints(SolutionSet solutionSet) {
		extremePoints = new double[numberOfObjectives][numberOfObjectives];
		for (int j = 0; j < numberOfObjectives; j++) {
			int index = -1;
			double min = Double.MAX_VALUE;

			for (int i = 0; i < solutionSet.size(); i++) {
				double asfValue = asfFunction(solutionSet.get(i), j);
				if (asfValue < min) {
					min = asfValue;
					index = i;
				}
			}

			for (int k = 0; k < numberOfObjectives; k++) {
				extremePoints[j][k] = solutionSet.get(index).getObjective(k);
			}
		}
	}

	private double asfFunction(Solution sol, int j) {
		double max = Double.MIN_VALUE;
		double epsilon = 1.0E-6;

		for (int i = 0; i < numberOfObjectives; i++) {

			double val = Math.abs(sol.getObjective(i) - idealPoint[i]);

			if (j != i) {
				val = val / epsilon;
			}

			if (val > max) {
				max = val;
			}
		}

		return max;
	}

	private void computeIntercepts() {

		intercepts = new double[numberOfObjectives];

		double[][] temp = new double[numberOfObjectives][numberOfObjectives];

		for (int i = 0; i < numberOfObjectives; i++) {
			for (int j = 0; j < numberOfObjectives; j++) {
				double val = extremePoints[i][j] - idealPoint[j];
				temp[i][j] = val;
			}
		}

		Matrix EX = new Matrix(temp);

		if (EX.rank() == EX.getRowDimension()) {
			double[] u = new double[numberOfObjectives];
			for (int j = 0; j < numberOfObjectives; j++) {
				u[j] = 1;
			}

			Matrix UM = new Matrix(u, numberOfObjectives);

			Matrix AL = EX.inverse().times(UM);

			int j = 0;
			for (j = 0; j < numberOfObjectives; j++) {

				double aj = 1.0 / AL.get(j, 0) + idealPoint[j];

				if ((aj > idealPoint[j]) && (!Double.isInfinite(aj))
						&& (!Double.isNaN(aj))) {
					intercepts[j] = aj;
				} else {
					break;
				}
			}
			if (j != numberOfObjectives) {
				for (int k = 0; k < numberOfObjectives; k++) {
					intercepts[k] = nadirPoint[k];
				}
			}

		} else {
			for (int k = 0; k < numberOfObjectives; k++) {
				intercepts[k] = nadirPoint[k];
			}
		}

	}

	private void normalizePopulation(SolutionSet solutionSet) {
		for (int i = 0; i < solutionSet.size(); i++) {
			Solution sol = solutionSet.get(i);

			for (int j = 0; j < numberOfObjectives; j++) {

				double val = (sol.getObjective(j) - idealPoint[j])
						/ (intercepts[j] - idealPoint[j]);
				sol.setTranslatedObjectives(val, j);
			}// for
		}

		for (int j = 0; j < numberOfObjectives; j++) {
			intercepts[j] = 1.0;
		}// for

	}

}
