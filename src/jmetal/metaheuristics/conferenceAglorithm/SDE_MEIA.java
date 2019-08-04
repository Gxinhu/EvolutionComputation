package jmetal.metaheuristics.conferenceAglorithm;

import jmetal.core.*;
import jmetal.metaheuristics.r2pso.util.ShiftedEuclideanDistanceAssigment;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.distanceToZmin;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Random;

public class SDE_MEIA extends Algorithm {
	private final double k = 0.06;
	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private QualityIndicator indicator;
	private int t;
	private int[][] neighborhood;
	private int populationSize;
	/**
	 * Stores the SolutionSet
	 */
	private NonDominatedSolutionList archive;
	private SolutionSet population;
	private SolutionSet tempPopulation;
	private SolutionSet clonePopulation;
	private SolutionSet lastPbestPopulation;
	/**
	 * Z vector (ideal point)
	 */

	double[][] angleMatrix;
	private boolean[] removed;
	private double[] minAngleArray;
	/**
	 * Lambda vectors
	 */
	private double[][] lambdaVectors;
	private double[] cosineLambda;
	private Random random = new Random();
	private int iteration;
	private Operator mutationOperator;
	private Operator crossoverOperator;
	private Operator cloneOperator;
	private double w;
	private double mu1, mu2, sigma1, sigma2;
	private int maxIterations;
	private double[] idealPoint;
	normalizationNSGAIII normalizationbynsga3;
	/* Calculate the Shifted Distance */

	private ShiftedEuclideanDistanceAssigment shiftedEuclideanDistanceAssigment;

	public SDE_MEIA(Problem problem) {
		super(problem);
		this.problem = problem;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		idealPoint = new double[problem.getNumberOfObjectives()];
		shiftedEuclideanDistanceAssigment = new ShiftedEuclideanDistanceAssigment(problem);
		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		iteration = 0;
		cosineLambda = new double[populationSize];
		archive = new NonDominatedSolutionList();
		population = new SolutionSet(populationSize);
		lastPbestPopulation = new SolutionSet(populationSize);
		tempPopulation = new SolutionSet(2 * populationSize);
		clonePopulation = new SolutionSet(2 * populationSize);
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		cloneOperator = operators_.get("clone");
		t = populationSize / 5;
		neighborhood = new int[populationSize][t];
		lambdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];
		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
		initPopulation();
		while (iteration < maxIterations) {
			cloneOffspringCreation();
			updatedAchieveByR2();
			++iteration;
		}
		return archive;
		//TODO is the population doesn't have a good diversity, so I think if population have a good diversity maybe can help the archive.
	}

	private void cloneOffspringCreation() throws JMException {
		// clone operation
		shiftedEuclideanDistanceAssigment.fitnessCompute(archive);
		archive.sort(new CrowdingComparator());
		int cloneSize = populationSize / 5;
		SolutionSet cPopulation = new SolutionSet(cloneSize);
		for (int k = 0; k < archive.size() && k < cloneSize; k++) {
			cPopulation.add(new Solution(archive.get(k)));
		}
		clonePopulation.clear();
		clonePopulation = (SolutionSet) cloneOperator.execute(cPopulation);
		tempPopulation.clear();
		for (int i = 0; i < clonePopulation.size(); i++) {
			Solution[] particle2 = new Solution[2];
			int ran;
			particle2[0] = clonePopulation.get(0);
			ran = PseudoRandom.randInt(0, clonePopulation.size() - 1);
			particle2[1] = clonePopulation.get(ran);
			Solution[] offSpring = (Solution[]) crossoverOperator.execute(particle2);
			mutationOperator.execute(offSpring[0]);
			problem.evaluate(offSpring[0]);
			tempPopulation.add(offSpring[0]);
		}

	}

	private void updatedAchieveByR2() {
		SolutionSet temp = new SolutionSet(2 * populationSize);
		for (int k = 0; k < tempPopulation.size(); k++) {
			archive.add(new Solution(tempPopulation.get(k)));
		}
		for (int i = 0; i < archive.size(); i++) {
			temp.add(new Solution(archive.get(i)));
		}
		if (temp.size() > populationSize) {
			archive.clear();
			normalizationbynsga3 = new normalizationNSGAIII(temp, problem.getNumberOfObjectives());
			normalizationbynsga3.execute();
			RealMatrix functionValueMatrix = new Array2DRowRealMatrix(temp.writeTransObjectivesToMatrix());
			int[] maxIndex = new int[problem.getNumberOfObjectives()];
			double[] maxValues = new double[problem.getNumberOfObjectives()];
			boolean[] removed = new boolean[temp.size()];
			boolean flag = true;
			while (flag) {
				flag = false;
				for (int j = 0; j < problem.getNumberOfObjectives(); j++) {
					maxIndex[j] = functionValueMatrix.getColumnVector(j).getMaxIndex();
					maxValues[j] = functionValueMatrix.getColumnVector(j).getMaxValue();
					if (maxValues[j] - 1e-6 > 1) {
						if (!removed[maxIndex[j]]) {
							removed[maxIndex[j]] = true;
							flag = true;
						}
					}
				}
			}
			for (int i = temp.size() - 1; i >= 0; i--) {
				if (removed[i]) {
					temp.remove(i);
				}
			}
			if (temp.size() > populationSize) {
				for (int i = 0; i < temp.size(); i++) {
					temp.get(i).setID(i);
				}
				normalizationbynsga3 = new normalizationNSGAIII(temp, problem.getNumberOfObjectives());
				normalizationbynsga3.execute();
				functionValueMatrix = new Array2DRowRealMatrix(temp.writeTransObjectivesToMatrix());
				calculateDistanceToZmin(temp, functionValueMatrix);
				initializeAngleMatrix(functionValueMatrix);
				eliminate(temp);
			} else {
				for (int i = 0; i < temp.size(); i++) {
					archive.add(new Solution(temp.get(i)));
				}
			}
		}
	}

	private void calculateDistanceToZmin(SolutionSet temp, RealMatrix functionValueMatrix) {
		double[] ideal = new double[problem.getNumberOfObjectives()];
		for (int i = 0; i < functionValueMatrix.getRowDimension(); i++) {
			RealVector solutionVector = functionValueMatrix.getRowVector(i);
			double tempDistance = solutionVector.getDistance(new ArrayRealVector(ideal));
			temp.get(i).setDistanceToZmin(tempDistance);
		}
	}

	/**
	 * @param union
	 * @param functionValueMatrix
	 */
	private void eliminate(SolutionSet union) {
		// Step: 1 Initialize angle matrix
		union.sort(new distanceToZmin());
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(union.writeTransObjectivesToMatrix());
		removed = new boolean[union.size()];
		minAngleArray = new double[union.size()];
		boolean[] isExtreme = new boolean[union.size()];
		boolean[] considered = new boolean[union.size()];
		int[] removedSolutions = new int[populationSize];
		int noOfRemoved = 0;
		// Add extreme solutions

		for (int k = 0; k < problem.getNumberOfObjectives(); k++) {
			double[] vector1 = new double[problem.getNumberOfObjectives()];
			vector1[k] = 1.0;
			RealVector vector = new ArrayRealVector(vector1);
			double minAngle = 1e+30;
			int minAngleId = -1;

			for (int i = 0; i < union.size(); i++) {

				if (!removed[i]) {
					RealVector solutionVector = functionValueMatrix.getRowVector(i);
					double angle = 1 - solutionVector.cosine(vector);
					if (angle < minAngle) {
						boolean flag = true;
//						for (int j = 0; j < problem.getNumberOfObjectives(); j++) {
						if (solutionVector.getMaxValue() - 1e-6 > 1) {
							flag = false;
						}
//						}
						if (flag) {
							minAngle = angle;
							minAngleId = i;
						}
					} // for i
				}
			}
			archive.add(new Solution(union.get(minAngleId)));
			removed[minAngleId] = true;
			isExtreme[minAngleId] = true;
			considered[minAngleId] = true;
		} // for k
		boolean removedFlag = false;

		// Find the minimum angle for each unadded solution
		for (int i = 0; i < union.size(); i++) {
			if (!removed[i]) {
				double minAng = 1.0e+30;
				int minAngID = -1;

				for (int j = 0; j < union.size(); j++) {

					if (j == i) {
						continue;
					}
					try {
						if (!removed[j]
								&& angleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
							minAng = angleMatrix[union.get(i).getID()][union.get(j).getID()];
							minAngID = j;
						}
					} catch (ArrayIndexOutOfBoundsException e) {
						e.printStackTrace();
					}

				} // for j
				union.get(i).setClusterID(minAngID);
				union.get(i).setAssociateDist(minAng);
				minAngleArray[i] = minAng;
			} // if
		} // for i

		int remain = populationSize - archive.size();
//		System.out.println("remain" + remain);
		//  Remove solutions
		int numberOfLoops = union.size() - archive.size() - remain;

		if (removedFlag) {
			numberOfLoops = (union.size() - problem.getNumberOfObjectives()) - archive.size() - remain;
		}
		for (int r = 0; r < numberOfLoops; r++) {

			double minValue = 1.0e+30;
			int minValueID = -1;

			for (int i = 0; i < union.size(); i++) {
				if (removed[i]) {
					continue;
				}
				if (minAngleArray[i] <= minValue) {
					minValue = minAngleArray[i];
					minValueID = i;
				}
			} // for i

			int associatedId = union.get(minValueID).getClusterID();
			int removedID = -1;
			if (union.get(minValueID).getDistanceToZmin() < union.get(associatedId).getDistanceToZmin()) {

				removedID = associatedId;

			} else if (union.get(minValueID).getDistanceToZmin() > union.get(associatedId).getDistanceToZmin()) {

				removedID = minValueID;

			} else {

				if (PseudoRandom.randDouble() < 0.5) {
					removedID = associatedId;
				} else {
					removedID = minValueID;
				}

//				}// IF

			} // if

			removed[removedID] = true;

			removedSolutions[noOfRemoved] = removedID;
			noOfRemoved++;

			considered[minValueID] = true;
			considered[associatedId] = true;

			// Update angles
			for (int i = 0; i < union.size(); i++) {
				if (removed[i]) {
					continue;
				}

				if (angleMatrix[union.get(i).getID()][union.get(removedID).getID()] == union.get(i).getAssociateDist()) {
					double minAng = 1.0e+30;
					int minAngId = -1;

					for (int j = 0; j < union.size(); j++) {

						if (j == i) {
							continue;
						}

						if (!removed[j] && angleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
							minAng = angleMatrix[union.get(i).getID()][union.get(j).getID()];
							minAngId = j;
						}

					} // for j

					union.get(i).setClusterID(minAngId);
					union.get(i).setAssociateDist(minAng);
					minAngleArray[i] = minAng;
				}
			}

		} // for r

		// Add remain solutions into pop

		for (int i = 0; i < union.size(); i++) {
			if (!removed[i] && !isExtreme[i]) {
				archive.add(new Solution(union.get(i)));
			}
		} // for i


	}

	private void initializeAngleMatrix(RealMatrix functionValueMatrix) {

		angleMatrix = new double[functionValueMatrix.getRowDimension()][functionValueMatrix.getRowDimension()];
		for (int i = 0; i < functionValueMatrix.getRowDimension(); i++) {
			for (int j = i; j < functionValueMatrix.getRowDimension(); j++) {
				if (i != j) {
					RealVector functionValueVector1 = functionValueMatrix.getRowVector(i);
					RealVector functionValueVector2 = functionValueMatrix.getRowVector(j);
					angleMatrix[i][j] = 1 - functionValueVector1.cosine(functionValueVector2);
					angleMatrix[j][i] = angleMatrix[i][j];
				} else {
					angleMatrix[i][j] = 0.0;
				}
			}
		}
	}

	public void initPopulation() throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem);
			problem.evaluate(newSolution);
			population.add(new Solution(newSolution));
			archive.add(new Solution(newSolution));
		}
	}

}