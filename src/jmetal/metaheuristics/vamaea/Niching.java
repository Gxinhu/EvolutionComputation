package jmetal.metaheuristics.vamaea;

import Jama.Matrix;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;

public class Niching {
	SolutionSet population;
	SolutionSet lastFront;
	SolutionSet mgPopulation;

	SolutionSet union;

	int obj;
	int remain;

	boolean normalization;

	double[] zideal;

	double[] zmax;

	private double[] zp_; // 	ideal point 
	private double[] nzp_; // nadir point

	double[][] extremePoints;

	double[] intercepts;

	int div;

	double alpha;

	int populationSize;

	public Niching(SolutionSet population, SolutionSet lastFront,
	               int remain, boolean normalization) {

		this.population = population;
		this.lastFront = lastFront;

		this.remain = remain;

		this.normalization = normalization;

		this.mgPopulation = population.union(lastFront);

		if (population.size() > 0) {
			this.obj = population.get(0).getNumberOfObjectives();
		} else {
			this.obj = lastFront.get(0).getNumberOfObjectives();
		}
	}

	public Niching(SolutionSet population, SolutionSet union,
	               boolean normalization, int div, int populationSize) {

		this.population = population;
		this.mgPopulation = union;
		this.normalization = normalization;
		this.div = div;
		this.populationSize = populationSize;

		if (population.size() > 0) {
			this.obj = population.get(0).getNumberOfObjectives();
		} else {
			this.obj = union.get(0).getNumberOfObjectives();
		}

	}


	public Niching(SolutionSet population, SolutionSet union,
	               boolean normalization, int populationSize) {

		this.population = population;
		this.mgPopulation = union;
		this.normalization = normalization;
		this.populationSize = populationSize;

		if (population.size() > 0) {
			this.obj = population.get(0).getNumberOfObjectives();
		} else {
			this.obj = union.get(0).getNumberOfObjectives();
		}

	}

	public Niching(SolutionSet population, SolutionSet union,
	               boolean normalization, int populationSize, double alpha) {

		this.population = population;
		this.mgPopulation = union;
		this.normalization = normalization;
		this.populationSize = populationSize;
		this.alpha = alpha;

		if (population.size() > 0) {
			this.obj = population.get(0).getNumberOfObjectives();
		} else {
			this.obj = union.get(0).getNumberOfObjectives();
		}

	}

	public Niching(SolutionSet population, SolutionSet union,
	               boolean normalization, int populationSize, double[] zp, double[] nzp) {

		this.population = population;
		this.mgPopulation = union;
		this.normalization = normalization;
		this.populationSize = populationSize;
		this.zp_ = zp;
		this.nzp_ = nzp;

		if (population.size() > 0) {
			this.obj = population.get(0).getNumberOfObjectives();
		} else {
			this.obj = union.get(0).getNumberOfObjectives();
		}

	}

	public void execute() {
		computeIdealPoint();

		if (normalization) {
			computeMaxPoint();

			/**
			 * Method 1 : The same as in NSGAIII
			 */
			computeExtremePoints();
			computeIntercepts();
			normalizePopulation();

			computeNorm(mgPopulation);
			rank();
		}
		assignment();
	}

	void computeIdealPoint() {
		zideal = new double[obj];

		for (int j = 0; j < obj; j++) {
			zideal[j] = 1.0e+30;

			for (int i = 0; i < mgPopulation.size(); i++) {
				if (mgPopulation.get(i).getObjective(j) < zideal[j]) {
					zideal[j] = mgPopulation.get(i).getObjective(j);
				}
			}
		}

	}

	void computeMaxPoint() {
		zmax = new double[obj];

		for (int j = 0; j < obj; j++) {
			zmax[j] = -1.0e+30;

			for (int i = 0; i < mgPopulation.size(); i++) {
				if (mgPopulation.get(i).getObjective(j) > zmax[j]) {
					zmax[j] = mgPopulation.get(i).getObjective(j);
				}
			}
		}
	}

	void computeExtremePoints() {
		extremePoints = new double[obj][obj];

		for (int j = 0; j < obj; j++) {
			int index = -1;
			double min = Double.MAX_VALUE;

			for (int i = 0; i < mgPopulation.size(); i++) {
				double asfValue = asfFunction(mgPopulation.get(i), j);
				if (asfValue < min) {
					min = asfValue;
					index = i;
				}
			}

			for (int k = 0; k < obj; k++) {
				extremePoints[j][k] = mgPopulation.get(index).getObjective(k);
			}
		}
	}

	void computeIntercepts() {

		intercepts = new double[obj];

		double[][] temp = new double[obj][obj];

		for (int i = 0; i < obj; i++) {
			for (int j = 0; j < obj; j++) {
				double val = extremePoints[i][j] - zideal[j];
				temp[i][j] = val;
			}
		}

		Matrix EX = new Matrix(temp);

		if (EX.rank() == EX.getRowDimension()) {
			double[] u = new double[obj];
			for (int j = 0; j < obj; j++) {
				u[j] = 1;
			}

			Matrix UM = new Matrix(u, obj);

			Matrix AL = EX.inverse().times(UM);

			int j = 0;
			for (j = 0; j < obj; j++) {

				double aj = 1.0 / AL.get(j, 0) + zideal[j];

				if ((aj > zideal[j]) && (!Double.isInfinite(aj)) && (!Double.isNaN(aj))) {
					intercepts[j] = aj;
				} else {
					break;
				}
			}
			if (j != obj) {
				for (int k = 0; k < obj; k++) {
					intercepts[k] = zmax[k];
				}
			}

		} else {
			for (int k = 0; k < obj; k++) {
				intercepts[k] = zmax[k];
			}
		}

	}

	void normalizePopulation() {
		for (int i = 0; i < mgPopulation.size(); i++) {
			Solution sol = mgPopulation.get(i);
			double fitness = 0.0;

			for (int j = 0; j < obj; j++) {

				double val = (sol.getObjective(j) - zideal[j])
						/ (intercepts[j] - zideal[j]);
				fitness = fitness + val;
				sol.setTranslatedObjectives(val, j);
			}// for

			sol.setFitness(fitness);
		}
	}


	void simpleNormalizePopulation() {
		for (int i = 0; i < mgPopulation.size(); i++) {
			Solution sol = mgPopulation.get(i);
			double fitness = 0.0;

			for (int j = 0; j < obj; j++) {
				double val = 0.0;

				val = (sol.getObjective(j) - zideal[j])
						/ (zmax[j] - zideal[j]);

				fitness = fitness + val;
				sol.setTranslatedObjectives(val, j);

			}// for 

			sol.setFitness(fitness);
		}
	}


	double asfFunction(Solution sol, int j) {
		double max = Double.MIN_VALUE;
		double epsilon = 1.0E-6;

		for (int i = 0; i < obj; i++) {

			double val = Math.abs(sol.getObjective(i) - zideal[i]);

			if (j != i) {
				val = val / epsilon;
			}

			if (val > max) {
				max = val;
			}
		}

		return max;
	}

	public void computeNorm(SolutionSet front) {
		// Calculate norm of each solution in front 		
		for (int i = 0; i < front.size(); i++) {
			Solution sol = front.get(i);

			double norm = 0.0;
			for (int j = 0; j < obj; j++) {
				norm += sol.getTranslatedObjectives(j) * sol.getTranslatedObjectives(j);
			}
			norm = Math.sqrt(norm);
			sol.setDistanceToZmin(norm); // This is the norm of the individual
		}
	}

	public void rank() {
		// Ranking the union
		Ranking ranking = new NondominatedRanking(mgPopulation);

		remain = populationSize;
		int index = 0;

		population.clear();

		// Obtain the next front
		lastFront = ranking.getSubfront(index);

		while ((remain > 0) && (remain >= lastFront.size())) {

			for (int k = 0; k < lastFront.size(); k++) {
				population.add(lastFront.get(k));
			} // for

			// Decrement remain
			remain = remain - lastFront.size();

			// Obtain the next front
			index++;
			if (remain > 0) {
				lastFront = ranking.getSubfront(index);
			} // if
		}

	}


	public void assignment() {

		if (remain > 0) {
			int n = lastFront.size();
			double[] angles = new double[n];
			int[] index = new int[n];
			boolean[] removed = new boolean[n];
			lastFront.sort(new FitnessComparator());

			/**
			 *  If the population is empty that is common in higher objective space
			 */
			if (population.size() == 0) {

				for (int o = 0; o < this.obj; o++) {
					double minAngle2Axis = 1.0e+30;
					int minAngle2AxisID = -1;

					for (int i = 0; i < n; i++) {
						if (removed[i] == false) {
							Solution solLastFront = lastFront.get(i);
							double angle = Math.acos(Math.abs(solLastFront.getTranslatedObjectives(o) / solLastFront.getDistanceToZmin()));

							if (angle < minAngle2Axis) {
								minAngle2Axis = angle;
								minAngle2AxisID = i;
							}
						}
					}// for 

					removed[minAngle2AxisID] = true;
					population.add(new Solution(lastFront.get(minAngle2AxisID)));
					remain--;
				} // for o 

				int k = 0;
				int t = 0;
				while (k < this.obj && remain > 0) {
					if (removed[t] == false) {
						population.add(new Solution(lastFront.get(t)));
						removed[t] = true;
						k++;
						remain--;
					}
					t++;
				} // Add better solutions in terms of the fitness value

			} // If population.size == 0

			/**
			 * Associate each solution in the last front with a solution in the population
			 */
			for (int i = 0; i < n; i++) {
				Solution solLastFront = lastFront.get(i);
				double minAng = 1.0e+30;
				int minAngID = -1;

				for (int j = 0; j < population.size(); j++) {
					Solution solPop = population.get(j);
					double angle = calAngle(solLastFront, solPop);
					if (angle < minAng) {
						minAng = angle;
						minAngID = j;
					}
				} // for j		
				angles[i] = minAng;
				index[i] = minAngID;
			} // for i		

			/**
			 * Niching procedure
			 */
			for (int r = 0; r < remain; r++) {
				/**
				 * Step 1: Find max and min angles 
				 */
				int maxAngleID = -1;
				double maxAngle = -1.0e+30;

				int minAglID = -1;
				double minAgl = 1.0e+30;

				for (int j = 0; j < n; j++) {
					// Find max angle 
					if (removed[j] == false && angles[j] > maxAngle) {
						maxAngle = angles[j];
						maxAngleID = j;
					}

					// Find min angle 
					if (removed[j] == false && angles[j] < minAgl) {
						minAgl = angles[j];
						minAglID = j;
					}
				} // for


				/**
				 * Step 2: Maximum-angle-first principle
				 *
				 */

				if (maxAngleID != -1) {    // Not all solutions in the last front have been added
					removed[maxAngleID] = true;
					population.add(new Solution(lastFront.get(maxAngleID)));

					// Update angles
					for (int i = 0; i < n; i++) { // For each solution in the last front
						if (removed[i] == false) {
							double angle = calAngle(lastFront.get(i), lastFront.get(maxAngleID));
							if (angle < angles[i]) {
								angles[i] = angle;
								index[i] = population.size() - 1;
							}
						}//if
					}//for i

				} else {
					break;
				}


				/**
				 * Step 2: Worse elimination principle
				 *
				 */

				if (minAglID != -1 && minAgl < Math.PI / 2 / (this.populationSize + 1)) {

					double norm1 = population.get(index[minAglID]).getFitness();
					double norm2 = lastFront.get(minAglID).getFitness();

					if (norm1 > norm2 && removed[minAglID] == false) {

						removed[minAglID] = true;
						population.replace(index[minAglID], new Solution(lastFront.get(minAglID)));

						// update angles
						for (int i = 0; i < n; i++) {
							if (removed[i] == false) {

								if (index[i] != index[minAglID]) {
									double angle = calAngle(lastFront.get(i), lastFront.get(minAglID));
									if (angle < angles[i]) {
										angles[i] = angle;
										index[i] = index[minAglID];
									}//if
								} else {
									double angle = calAngle(lastFront.get(i), lastFront.get(minAglID));
									if (angle < angles[i]) {
										angles[i] = angle;
										index[i] = index[minAglID];
									} else {
										// find the minimum angle
										double minAng = 1.0e+30;
										int minAngID = -1;

										for (int j = 0; j < population.size(); j++) {
											Solution solPop = population.get(j);
											angle = calAngle(lastFront.get(i), solPop);
											if (angle < minAng) {
												minAng = angle;
												minAngID = j;
											}
										} // for j		
										angles[i] = minAng;
										index[i] = minAngID;
									}

								}
							}//if
						}//for i
					}
				}    // if minAglID

				/*----------------------------------------Line ---------------------------------------------*/
			}// for r 
		}//if remain > 0
	}


	public double calAngle(Solution s1, Solution s2) {
		double angle = 0.0;
		double norm1 = 0.0;
		double norm2 = 0.0;

		double innerProduct = 0.0;

		for (int i = 0; i < obj; i++) {
			innerProduct += s1.getTranslatedObjectives(i) * s2.getTranslatedObjectives(i);
		}

		norm1 = s1.getDistanceToZmin();
		norm2 = s2.getDistanceToZmin();

		angle = Math.acos(Math.abs(innerProduct / (norm1 * norm2)));

		return angle;
	}

	public double calAngle(Solution s1, double[] vector) {
		double angle = 0.0;
		double norm1 = 0.0;
		double norm2 = 1.0;

		double innerProduct = 0.0;

		for (int i = 0; i < obj; i++) {
			innerProduct += s1.getTranslatedObjectives(i) * vector[i];
		}

		norm1 = s1.getDistanceToZmin();

		angle = Math.acos(Math.abs(innerProduct / (norm1 * norm2)));

		return angle;
	}

	static void QuickSort(double[] array, int[] idx, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int tempIdx = idx[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
					int tempIndex = idx[j];
					idx[j] = idx[i];
					idx[i] = tempIndex;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			idx[to] = idx[i + 1];
			idx[i + 1] = tempIdx;
			QuickSort(array, idx, from, i);
			QuickSort(array, idx, i + 1, to);
		}
	}
}
