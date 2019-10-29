package jmetal.metaheuristics.Vapso;

import jmetal.core.*;
import jmetal.metaheuristics.Vapso.util.ShiftedEuclideanDistanceAssigment;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.distanceToZmin;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.Random;

public class VagPsodeletezmin extends Algorithm {
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

	public VagPsodeletezmin(Problem problem) {
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
		t = 20;
		neighborhood = new int[populationSize][t];
		lambdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];
		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
		initNeighborhood();
		initPopulation();
		orderPopulation(population);
		while (iteration < maxIterations) {
			cloneOffspringCreation();
			updatedAchieveByR2();
			++iteration;
			offspringcreationbypso();
			updatedAchieveByR2();
			++iteration;
		}
		return archive;
		//TODO is the population doesn't have a good diversity, so I think if population have a good diversity maybe can help the archive.
	}

	private void orderPopulation(SolutionSet pop) {
		population = new SolutionSet(populationSize);

		double[][] fitnesses = new double[this.populationSize][this.populationSize];
		double[][] objective = pop.writeObjectivesToMatrix();
		for (int i = 0; i < this.populationSize; i++) {
			for (int j = 0; j < this.populationSize; j++) {
				fitnesses[i][j] = this.pbi(objective[i],
						this.lambdaVectors[j], 5);
			}
		}
		for (int i = 0; i < this.populationSize; i++) {
			double minFit = Double.MAX_VALUE;
			int particleIndex = -1;
			for (int j = 0; j < this.populationSize; j++) {
				if (fitnesses[j][i] < minFit) {
					minFit = fitnesses[j][i];
					particleIndex = j;
				}
			}
			this.population.add(pop.get(particleIndex));
			for (int n = 0; n < this.populationSize; n++) {
				fitnesses[particleIndex][n] = Double.MAX_VALUE;
			}
			fitnesses[particleIndex][i] = Double.MAX_VALUE;
			//this.leaders_.add(pop.get(particleIndex));
			this.archive.add(pop.get(particleIndex));
		}

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
			int removeNumber = 0;
			while (flag & (temp.size() - removeNumber >= populationSize)) {
				flag = false;
				for (int j = 0; j < problem.getNumberOfObjectives(); j++) {
					maxIndex[j] = functionValueMatrix.getColumnVector(j).getMaxIndex();
					maxValues[j] = functionValueMatrix.getColumnVector(j).getMaxValue();
					if (maxValues[j] - 1e-6 > 1) {
						removed[maxIndex[j]] = true;
						functionValueMatrix.setEntry(maxIndex[j], j, -1);
						flag = true;
						removeNumber += 1;
					}
				}
			}

			for (int i = temp.size() - 1; i >= 0; i--) {
				if (removed[i]) {
					temp.remove(i);
				}
			}
			for (int j = 0; j < problem.getNumberOfObjectives(); j++) {
				maxIndex[j] = functionValueMatrix.getColumnVector(j).getMaxIndex();
				maxValues[j] = functionValueMatrix.getColumnVector(j).getMaxValue();
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
	 */
	private void eliminate(SolutionSet union) {
		// Step: 1 Initialize angle matrix
		union.sort(new distanceToZmin());
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(union.writeTransObjectivesToMatrix());
		removed = new boolean[union.size()];
		minAngleArray = new double[union.size()];
		boolean[] isExtreme = new boolean[union.size()];
		boolean[] considered = new boolean[union.size()];
		int[] removedSolutions = new int[(int) (populationSize * 1.5)];
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
						minAngle = angle;
						minAngleId = i;
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

		int remain = (int) (0.9 * populationSize) - archive.size();
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
		SolutionSet removePopulation = new SolutionSet(2 * populationSize);
		for (int i = 0; i < union.size(); i++) {
			if (!removed[i] && !isExtreme[i]) {
				archive.add(new Solution(union.get(i)));
			} else if (removed[i] && !isExtreme[i]) {
				removePopulation.add(new Solution(union.get(i)));
			}
		} // for i
		for (int i = 0; i < removePopulation.size(); i++) {
			double minAng = 1.0e+30;
			int minAngID = -1;
			for (int j = problem.getNumberOfObjectives(); j < archive.size(); j++) {
				if (angleMatrix[removePopulation.get(i).getID()][archive.get(j).getID()] < minAng) {
					minAng = angleMatrix[removePopulation.get(i).getID()][archive.get(j).getID()];
					minAngID = j;
				}


			} // for j
			removePopulation.get(i).setClusterID(minAngID);
			removePopulation.get(i).setAssociateDist(minAng);
			minAngleArray[i] = minAng;
		} // if
		remain = populationSize - archive.size();
		boolean[] deletemark = new boolean[removePopulation.size()];
		int discard = removePopulation.size() - remain;
		int discardss = 0;
		for (int i = 0; i < remain; i++) {
			int maxId = new ArrayRealVector(minAngleArray).getMaxIndex();
			int maxIndex = removePopulation.get(maxId).getID();
			archive.add(new Solution(removePopulation.get(maxId)));
			minAngleArray[maxId] = Double.NaN;
			deletemark[maxId] = true;
			for (int j = 0; j < removePopulation.size(); j++) {
				if (!deletemark[j]) {
					double minAngle = angleMatrix[maxIndex][removePopulation.get(j).getID()];
					if (minAngle < minAngleArray[j]) {
						minAngleArray[j] = minAngle;
						removePopulation.get(j).setClusterID(archive.size() - 1);
					}
				}
			}
			int minID = new ArrayRealVector(minAngleArray).getMinIndex();
			int minIndex = removePopulation.get(minID).getID();
			int associateID = removePopulation.get(minID).getClusterID();
			try {
				if (archive.get(associateID).getDistanceToZmin() > removePopulation.get(minID).getDistanceToZmin()) {
					discardss++;
					minAngleArray[minID] = Double.NaN;
					deletemark[minID] = true;
					archive.replace(associateID, new Solution(removePopulation.get(minID)));
				}
			} catch (IndexOutOfBoundsException e) {
				e.printStackTrace();
			}
			for (int j = 0; j < removePopulation.size(); j++) {
				if (!deletemark[j]) {
					double minAngle = angleMatrix[minIndex][removePopulation.get(j).getID()];
					if (minAngle < minAngleArray[j]) {
						minAngleArray[j] = minAngle;
						removePopulation.get(j).setClusterID(associateID);
					}
				}
			}
			if (discardss >= discard) {
				i = remain;
			}
		}
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

	/**
	 * Modify the Coefficient used in Pso
	 *
	 * @param index
	 * @param angle
	 */
	private void calculateCoefficientValuesModify(int index, double[] angle) {
		//calculate Vc
		double vc;
		double distance = 0;
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			distance += Math.pow(idealPoint[i] - population.get(index).getObjective(i), 2);
		}
		distance = Math.sqrt(distance);
		vc = distance;
		double p1;
//		double theta=0.2-Math.acos(angle[index]);
		double theta = (Math.PI / 2 - angle[index]) / (Math.PI / 2) - 0.5;
//		theta=-0.5;
		p1 = theta;
		//calculate F
		double f;
		f = 0.25;
		double a = Math.sqrt(f);
		double m1 = (a + 1.0) * (a + 1.0) * (a * a + 3.0 * a + 1.0);
		double m2 = (a + 1.0) * (a + 1.0) * (2.0 * a * a + 3.0 * a + 2.0);
		w = (m1 * vc + m2 * p1 * vc + p1 - 1.0) / (m2 * vc + m1 * p1 * vc - p1 + 1.0);
		double c = 2.0 * (1 - p1) * (w + 1.0) / (a + 1.0);
		mu1 = c / 2.0;
		mu2 = a * c / 2.0;
		sigma1 = c / Math.sqrt(12.0);
		sigma2 = a * c / Math.sqrt(12.0);
	}

	private void offspringcreationbypso() throws JMException {
		double[][] functionValueMatrix;
		normalizationbynsga3 = new normalizationNSGAIII(archive, problem.getNumberOfObjectives());
		normalizationbynsga3.execute();
		if (archive.size() != 1) {
			functionValueMatrix = archive.writeTransObjectivesToMatrix();
		} else {
			functionValueMatrix = archive.writeObjectivesToMatrix();
		}

		normalizationbynsga3 = new normalizationNSGAIII(population, problem.getNumberOfObjectives());
		normalizationbynsga3.execute();
		int[] pbestIndex = new int[populationSize];
		int best_ind;
		double minFit;
		RealMatrix populationValuesMatrix = new Array2DRowRealMatrix(population.writeTransObjectivesToMatrix());
		int[] thetaId = new int[populationSize];
		double[] angle = new double[populationSize];
		double[] theta = new double[populationSize];
		RealMatrix lambdaMatrix = new Array2DRowRealMatrix(lambdaVectors);
//		for (int i = 0; i < this.populationSize; i++) {
//			RealVector populationVector = populationValuesMatrix.getRowVector(i);
//			if (populationVector.getNorm() == 0) {
//				thetaId[i] = PseudoRandom.randInt(0, population.size() - 1);
//				angle[i] = 0;
//			} else {
//				double minAngle = Double.NEGATIVE_INFINITY;
//				int minId = -1;
//				for (int j = 0; j < populationSize; j++) {
//					RealVector lambdaVector = lambdaMatrix.getRowVector(j);
//					double tempAngle = lambdaVector.cosine(populationVector);
//					if (tempAngle > minAngle) {
//						minAngle = tempAngle;
//						minId = j;
//					}
//				}
//				if (minId == -1) {
//					minId = PseudoRandom.randInt(0, population.size() - 1);
//					minAngle = 0;
//				}
//				thetaId[i] = minId;
//				angle[i] = minAngle;
//			}
//		}
		for (int i = 0; i < populationSize; i++) {
			RealVector lambdaVector = lambdaMatrix.getRowVector(i);
			RealVector solutionVector = populationValuesMatrix.getRowVector(i);
			if (solutionVector.getNorm() == 0) {
				thetaId[i] = i;
				angle[i] = 0;
			} else {
				double tempangle = lambdaVector.cosine(solutionVector);
				thetaId[i] = i;
				angle[i] = Math.acos(tempangle);
			}
		}


		for (int i = 0; i < populationSize; i++) {
			theta[i] = k * problem.getNumberOfObjectives() * (Math.toDegrees(angle[i]));
		}
		double fitnesse;
		for (int i = 0; i < this.populationSize; i++) {
			best_ind = -1;
			minFit = Double.MAX_VALUE;

			for (int j = 0; j < archive.size(); j++) {
				fitnesse = this.pbi(functionValueMatrix[j], lambdaVectors[thetaId[i]], theta[i]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					best_ind = j;
				}
			}
			if (best_ind == -1) {
				best_ind = 0;
			}
			pbestIndex[i] = best_ind;
		}
		double[] angleArchive = new double[populationSize];
		for (int i = 0; i < populationSize; i++) {
			RealVector vector1 = new ArrayRealVector(lambdaVectors[thetaId[i]]);
			RealVector vector2 = new ArrayRealVector(functionValueMatrix[pbestIndex[i]]);
			if (Double.isNaN(vector2.getNorm())) {
				System.out.println("Hello");
			}
			angleArchive[i] = Math.acos(vector1.cosine(vector2));
		}
		updatePopulationPso(pbestIndex, angleArchive, thetaId);
	}

	private void updatePopulationPso(int[] pbestIndex, double[] angle, int[] thetaId) throws JMException {
		int ran;
		Variable[] pbest, gbest;
		SolutionSet temps = population.union(archive);
		RealMatrix tempMatrix = new Array2DRowRealMatrix(temps.writeObjectivesToMatrix());
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = tempMatrix.getColumnVector(i).getMinValue();
		}
		for (int i = 0; i < populationSize; i++) {
			calculateCoefficientValuesModify(i, angle);
			Variable[] paritcle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			pbest = archive.get(pbestIndex[i]).getDecisionVariables();
			ran = PseudoRandom.randInt(0, t - 1);
			gbest = archive.get(pbestIndex[neighborhood[thetaId[i]][ran]]).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				double theta1 = random.nextGaussian() * sigma1 + mu1;
				double theta2 = random.nextGaussian() * sigma2 + mu2;
				double temp = (w * velocity[j]) + theta1 * (pbest[j].getValue() - paritcle[j].getValue())
						+ theta2 * (gbest[j].getValue() - paritcle[j].getValue());
				population.get(i).setSpeed(j, temp);
			}
		}
		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle = ;
			for (int var = 0; var < population.get(n).getDecisionVariables().length; var++) {
				population.get(n).getDecisionVariables()[var].setValue(population.get(n).getDecisionVariables()[var].getValue() + population.get(n).getSpeed()[var]);
				if (population.get(n).getDecisionVariables()[var].getValue() < problem.getLowerLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getLowerLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
				if (population.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
			}
			problem.evaluate(population.get(n));
		}
		tempPopulation.clear();
		for (int i = 0; i < population.size(); i++) {
			tempPopulation.add(new Solution(population.get(i)));
		}
	}

	public void initNeighborhood() {
		RealMatrix lambdaMatrix = new Array2DRowRealMatrix(lambdaVectors);
		lambdaMatrix = lambdaMatrix.multiply(lambdaMatrix.transpose());
		for (int k = 0; k < populationSize; k++) {
			double[] arrays = lambdaMatrix.getRowVector(k).toArray();
			ArrayList<Integer> index = new ArrayList<>(arrays.length);
			for (int i = 0; i < arrays.length; i++) {
				index.add(i);
			}
			index.sort((o1, o2) -> Double.compare(arrays[o2], arrays[o1]));
			for (int j = 0; j < populationSize; j++) {
				if (j < t) {
					neighborhood[k][j] = index.get(j);
				}
			}
			cosineLambda[k] = Math.acos(arrays[index.get(1)]);
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

	private double pbi(double[] indiv, double[] lambda, double theta) {
		int i;
		double d1, d2, nl;
		double fin;

		d1 = d2 = nl = 0.0;
		for (i = 0; i < problem.getNumberOfObjectives(); i++) {
			d1 += (indiv[i]) * lambda[i];
			nl += Math.pow(lambda[i], 2.0);
		}
		d1 = Math.abs(d1) / Math.sqrt(nl);
		if (nl == 0.0) {
			System.out
					.println("ERROR: dived by zero(bad weihgted vector)\n");
			System.exit(0);
		}
		for (i = 0; i < problem.getNumberOfObjectives(); i++) {
			d2 += Math.pow((indiv[i])
					- (d1 * lambda[i]), 2.0);
		}
		d2 = Math.sqrt(d2);
		fin = (d1 + theta * d2);
		return fin;
	}
}