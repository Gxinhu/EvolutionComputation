package jmetal.metaheuristics.ragpso;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
import jmetal.util.PseudoRandom;
import jmetal.util.createWeight;
import jmetal.util.deepcopy.deepCopy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class ragpso extends Algorithm {
	final double k = 0.06;
	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private QualityIndicator indicator;
	private int t;
	private int[][] neighborhood;
	private int populationSize;
	/**
	 * Stores the SlutionSett
	 */
	private SolutionSet archive;
	private SolutionSet population;
	private SolutionSet tempPopulation;
	/**
	 * Z vector (ideal point)
	 */
	private double[] idealPoint;
	private double[] nadirPoint;
	/**
	 * Lambda vectors
	 */
	private double[][] lambdaVectors;
	private double[][] lambdaVectors0;
	private double[] cosineLambda;
	private Random random = new Random();
	private int runtime;
	private int iteration, iteration1, iteration2;
	private double w;
	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;
	private int maxIterations;

	ragpso(Problem problem, QualityIndicator indicator, int i) {
		super(problem);
		this.problem = problem;
		this.indicator = indicator;
		this.runtime = i;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		vMax = 25.0;
		vMin = 5.0;
		pMax = 0.8;
		pMin = 0.1;
		fMax = 25.0;
		fMin = 0.25;
		iteration1 = maxIterations / 5;
		iteration2 = 4 * maxIterations / 5;
		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		iteration = 0;
		cosineLambda = new double[populationSize];
		archive = new SolutionSet(populationSize);
		population = new SolutionSet(populationSize);
		tempPopulation = new SolutionSet(2 * populationSize);
		t = populationSize / 10;
		neighborhood = new int[populationSize][t];
		idealPoint = new double[problem.getNumberOfObjectives()];
		nadirPoint = new double[problem.getNumberOfObjectives()];
		lambdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];
		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
		lambdaVectors0 = deepCopy.deepCopysDouble2d(lambdaVectors);
		initNeighborhood();
		initPopulation();
		while (iteration < maxIterations) {
			calculateCoefficientValues();
			offspringCreationbypso();
			referenceselectpso();
			weightVectorAdaption();
			++iteration;
		}
		return archive;
	}

	private void calculateCoefficientValues() {
		//calculate Vc
		double vc;
		if (iteration < iteration1) {
			vc = vMax;
		} else if (iteration > iteration1 & iteration < iteration2) {
			vc = (iteration - iteration1) * (vMin - vMax) / (iteration2 - iteration1) + vMax;
		} else {
			vc = vMin;
		}
		//calculate p1
		double p1;
		if (iteration < iteration1) {
			p1 = pMin;
		} else if (iteration > iteration1 & iteration < (iteration2 - iteration1) / 2.0) {
			p1 = (iteration - iteration1) * (pMax - pMin) / ((double) (iteration2 - iteration1) / 2.0 - iteration1) + pMin;
		} else if (iteration < iteration2 & iteration > (iteration2 - iteration1) / 2.0) {
			p1 = (iteration - iteration1) * (pMin - pMax) / (iteration2 - (double) (iteration2 - iteration1) / 2.0) + pMax;
		} else {
			p1 = pMax;
		}
		//calculate F
		double f;
		if (iteration < iteration1) {
			f = fMin;
		} else if (iteration > iteration1 & iteration < iteration2) {
			f = 1;
		} else {
			f = fMax;
		}
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

	private void referenceselectpso() {
		SolutionSet tempAchieve = new SolutionSet(populationSize);
		//Objective value Translation
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			double minValue1 = functionValueMatrix.getColumnVector(i).getMinValue();
			double minValue2 = functionValueMatrix.getColumnVector(i).getMinValue();
			if (minValue1 < minValue2) {
				idealPoint[i] = minValue1;
			} else {
				idealPoint[i] = minValue2;
			}
		}
		for (int i = 0; i < archive.size(); i++) {
			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
		}
		//Population Partition
		//Calculate the cosine between referenceVector and X
		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
		double[][] cosine = new double[archive.size()][populationSize];
		for (int i = 0; i < archive.size(); i++) {
			for (int j = 0; j < populationSize; j++) {
				RealVector functionValueVector = functionValueMatrix.getRowVector(i);
				RealVector referenceVector = referenceMatrix.getRowVector(j);
				cosine[i][j] = functionValueVector.cosine(referenceVector);
			}
		}
		//Partition
		List[] subRegion = new List[populationSize];
		for (int i = 0; i < populationSize; i++) {
			subRegion[i] = new ArrayList();
		}
		RealMatrix cosineMatrix = new Array2DRowRealMatrix(cosine);
		int maxIndex;
		for (int i = 0; i < archive.size(); i++) {
			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
			subRegion[maxIndex].add(i);
		}
		//Angle-Penalized Distance(APD) Calculation
		for (int i = 0; i < populationSize; i++) {
			if (subRegion[i].size() != 0) {
				tempAchieve.add(chooseSolution(population.get(i), subRegion[i], i));
			} else {
				tempAchieve.add(new Solution(population.get(i)));
			}
		}
		archive.clear();
		for (int i = 0; i < tempAchieve.size(); i++) {
			archive.add(new Solution(tempAchieve.get(i)));
		}
	}

	private Solution chooseSolution(Solution solution, List list, int index) {
		NonDominatedSolutionList tempSolutions = new NonDominatedSolutionList(new jmetal.util.comparators.DominanceComparator());
		tempSolutions.add(solution);
		for (Object o : list) {
			tempSolutions.add(archive.get((Integer) o));
		}
		if (tempSolutions.size() == 1) {
			return tempSolutions.get(0);
		} else {
			RealMatrix functionValueMatrix = new Array2DRowRealMatrix(tempSolutions.writeObjectivesToMatrix());
			RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
			for (int i = 0; i < tempSolutions.size(); i++) {
				functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
			}
			RealVector functionValueVector = functionValueMatrix.getRowVector(0);
			RealVector referenceVector = referenceMatrix.getRowVector(index);
			double cosine = functionValueVector.cosine(referenceVector);
			double runIteration = Math.pow((double) iteration / maxIterations, 2);
			double theta = Math.acos(cosine) / cosineLambda[index];
			double norm = functionValueVector.getNorm();
			double minapd1 = (1 + problem.getNumberOfObjectives() * runIteration *
					theta) * norm;
			int minindex = 0;
			for (int j = 1; j < tempSolutions.size(); j++) {
				functionValueVector = functionValueMatrix.getRowVector(j);
				cosine = functionValueVector.cosine(referenceVector);
				theta = Math.acos(cosine) / cosineLambda[index];
				norm = functionValueVector.getNorm();
				double minpbi2 = (1 + problem.getNumberOfObjectives() * runIteration *
						theta) * norm;
				if (minpbi2 < minapd1) {
					minapd1 = minpbi2;
					minindex = j;
				}
			}
			return tempSolutions.get(minindex);
		}
	}

	private void offspringCreationbypso() throws JMException {
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
		RealMatrix pfunctionValuesMatrix = new Array2DRowRealMatrix(population.writeObjectivesToMatrix());
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = (functionValueMatrix.getColumnVector(i).getMinValue() < pfunctionValuesMatrix.getColumnVector(i).getMinValue())
					? functionValueMatrix.getColumnVector(i).getMinValue() : pfunctionValuesMatrix.getColumnVector(i).getMinValue();

		}
		for (int i = 0; i < archive.size(); i++) {
			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
			pfunctionValuesMatrix.setRowVector(i, pfunctionValuesMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
		}
		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
		double[][] cosine = new double[archive.size()][populationSize];
		for (int i = 0; i < archive.size(); i++) {
			for (int j = 0; j < populationSize; j++) {
				RealVector functionValueVector = functionValueMatrix.getRowVector(i);
				RealVector referenceVector = referenceMatrix.getRowVector(j);
				cosine[i][j] = functionValueVector.cosine(referenceVector);
			}
		}
		//Partition
		List[] subRegion = new List[populationSize];
		for (int i = 0; i < populationSize; i++) {
			subRegion[i] = new ArrayList();
		}
		RealMatrix cosineMatrix = new Array2DRowRealMatrix(cosine);
		int maxIndex;
		double[] angle = new double[populationSize];
		for (int i = 0; i < archive.size(); i++) {
			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
			angle[i] = cosineMatrix.getRowVector(i).getMaxValue();
			subRegion[maxIndex].add(i);
		}
		int[] pBestIndex = findpbest(subRegion, cosine, functionValueMatrix);
		updatePopulationPso(pBestIndex);

	}

	private int[] findpbest(List[] subRegion, double[][] cosine, RealMatrix functionValueMatrix) {
		int[] pbestIndex = new int[populationSize];
		for (int i = 0; i < populationSize; i++) {
			//1)当前权值向量为空那个么在周围邻居里面寻找一个作为pbset 当然要是周围寻找的还是为空的话就通过邻居的远近距离找个最近的，或者在所有邻居里面找然后通过自适应pbi来寻找pbest
			if (subRegion[i].size() == 1) {
				pbestIndex[i] = (int) subRegion[i].get(0);
			}
			//2）当前去去权值向量里面的subregion里面就只有一个个体那么就直接选择这个个体作为pbest
			else if (subRegion[i].size() == 0) {
				pbestIndex[i] = -1;
			}
			//3）不止一个那么就用自适应权值来比较哪一个好
			else {
				pbestIndex[i] = chooseSolutionIndex(subRegion[i], i, cosine, functionValueMatrix);
			}
		}
		for (int i = 0; i < populationSize; i++) {
			if (pbestIndex[i] == -1) {
				for (int j = 0; j < neighborhood[i].length; j++) {
					if (pbestIndex[neighborhood[i][j]] != -1) {
						pbestIndex[i] = pbestIndex[neighborhood[i][j]];
						break;
					}
				}
				if (pbestIndex[i] == -1) {
					pbestIndex[i] = PseudoRandom.randInt(0, archive.size());
				}
			}
		}
		return pbestIndex;
	}

	private int chooseSolutionIndex(List list, int index, double[][] cosine, RealMatrix functionValueMatrix) {
		NonDominatedSolutionList tempSolutions = new NonDominatedSolutionList(new jmetal.util.comparators.DominanceComparator());
		for (Object o : list) {
			tempSolutions.add(archive.get((Integer) o));
		}
		if (tempSolutions.size() == 1) {
			return (int) list.get(0);
		} else {
			double runIteration = Math.pow((double) iteration / maxIterations, 2);
			double theta = Math.acos(cosine[(int) list.get(0)][index]) / cosineLambda[index];
			double norm = functionValueMatrix.getRowVector((int) list.get(0)).getNorm();
			double minapd1 = (1 + problem.getNumberOfObjectives() * runIteration *
					theta) * norm;
			int minindex = 0;
			for (int j = 1; j < tempSolutions.size(); j++) {
				theta = Math.acos(cosine[(int) list.get(j)][index]) / cosineLambda[index];
				norm = functionValueMatrix.getRowVector((int) list.get(j)).getNorm();
				double minapd2 = (1 + problem.getNumberOfObjectives() * runIteration *
						theta) * norm;
				if (minapd2 < minapd1) {
					minapd1 = minapd2;
					minindex = j;
				}
			}
			return (int) list.get(minindex);
		}
	}

	private void updatePopulationPso(int[] pBestIndex) throws JMException {
		int ran;
		Variable[] pBest, gBest;
		for (int i = 0; i < populationSize; i++) {
			Variable[] particle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			pBest = archive.get(pBestIndex[i]).getDecisionVariables();
			ran = PseudoRandom.randInt(0, archive.size() - 1);
			gBest = archive.get(ran).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				double theta1 = random.nextGaussian() * sigma1 + mu1;
				double theta2 = random.nextGaussian() * sigma2 + mu2;
				double temp = (w * velocity[j]) + theta1 * (pBest[j].getValue() - particle[j].getValue())
						+ theta2 * (gBest[j].getValue() - particle[j].getValue());
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
	}

	private void weightVectorAdaption() {
		if (iteration % Math.ceil(maxIterations * 0.1) == 0) {
			RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
			RealVector subtract = new ArrayRealVector(problem.getNumberOfObjectives());
			for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
				idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
				nadirPoint[i] = functionValueMatrix.getColumnVector(i).getMaxValue();
				subtract.setEntry(i, nadirPoint[i] - idealPoint[i]);
			}
			RealMatrix lambdaMatrix = new Array2DRowRealMatrix(lambdaVectors0);
			for (int i = 0; i < populationSize; i++) {
				lambdaMatrix.setRowVector(i, lambdaMatrix.getRowVector(i).ebeMultiply(subtract));
				lambdaMatrix.setRowVector(i, lambdaMatrix.getRowVector(i).mapDivide(lambdaMatrix.getRowVector(i).getNorm()));
			}
			this.lambdaVectors = lambdaMatrix.getData();
			this.initNeighborhood();
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
			if (this.problem.getNumberOfConstraints() != 0) {
				problem.evaluateConstraints(newSolution);
			}
			// evaluations++;
			population.add(new Solution(newSolution));
			archive.add(new Solution(newSolution));
		}
	}

}