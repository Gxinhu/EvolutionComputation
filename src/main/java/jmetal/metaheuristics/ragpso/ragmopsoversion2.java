package jmetal.metaheuristics.ragpso;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;
import jmetal.util.deepcopy.deepCopy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.List;

public class ragmopsoversion2 extends Algorithm {
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
	private int runtime;
	private int iteration;
	private Operator mutationOperator;
	private Operator crossoverOperator;

	private int maxIterations;

	ragmopsoversion2(Problem problem, QualityIndicator indicator, int i) {
		super(problem);
		this.problem = problem;
		this.indicator = indicator;
		this.runtime = i;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		iteration = 0;
		cosineLambda = new double[populationSize];
		archive = new SolutionSet(populationSize);
		population = new SolutionSet(populationSize);
		tempPopulation = new SolutionSet(2 * populationSize);
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		t = populationSize / 5;
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
			offspringCreation();
			referenceSelect();
			weightVectorAdaption();
			++iteration;
			offspringCreationbyPSO();
			referenceSelectPSO();
			weightVectorAdaption();
			++iteration;
		}
		return archive;
	}

	private void referenceSelectPSO() {
		SolutionSet tempAchive = new SolutionSet(populationSize);
		//Objective value Translation
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
		RealMatrix pfuncitionValuesMatrix = new Array2DRowRealMatrix(population.writeObjectivesToMatrix());
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
		for (int i = 0; i < population.size(); i++) {
			pfuncitionValuesMatrix.setRowVector(i, pfuncitionValuesMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
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
				double degree = pfuncitionValuesMatrix.getRowVector(i).cosine(new ArrayRealVector(lambdaVectors[i]));
				double theta = k * problem.getNumberOfObjectives() * cosineLambda[i] * Math.acos(degree);
				tempAchive.add(chooseSolution(population.get(i), subRegion[i], theta, i));
			} else {
				tempAchive.add(new Solution(population.get(i)));
			}
		}
		archive = tempAchive;
	}

	private Solution chooseSolution(Solution solution, List list, double theta, int index) {
		NonDominatedSolutionList tempSolutionset = new NonDominatedSolutionList(new jmetal.util.comparators.DominanceComparator());
		tempSolutionset.add(solution);
		for (int i = 0; i < list.size(); i++) {
			tempSolutionset.add(archive.get((Integer) list.get(i)));
		}
		if (tempSolutionset.size() == 1) {
			return tempSolutionset.get(0);
		} else {
			double minpbi1 = pbi(tempSolutionset.get(0), lambdaVectors[index], theta);
			int minindex = 0;
			for (int j = 1; j < tempSolutionset.size(); j++) {
				double minpbi2 = pbi(tempSolutionset.get(j), lambdaVectors[index], theta);
				if (minpbi2 < minpbi1) {
					minpbi1 = minpbi2;
					minindex = j;
				}
			}
			return tempSolutionset.get(minindex);
		}
	}

	private void offspringCreationbyPSO() throws JMException {
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
		int[] pbestIndex = findPbest(subRegion, pfunctionValuesMatrix);
		updatePopulationPso(pbestIndex);

	}

	private int[] findPbest(List[] subRegion, RealMatrix pfunctionValuesMatrix) {
		int[] pbestIndex = new int[populationSize];
		for (int i = 0; i < populationSize; i++) {
			//1)当前权值向量为空那个么在周围邻居里面寻找一个作为pbset 当然要是周围寻找的还是为空的话就通过邻居的远近距离找个最近的，或者在所有邻居里面找然后通过自适应pbi来寻找pbest
			if (subRegion[i].size() == 1) {
				pbestIndex[i] = (int) subRegion[i].get(0);
			}
			//2）当前去去权值向量里面的subregion里面就只有一个个体那么就直接选择这个个体作为pbest
			else if (subRegion[i].size() == 0) {
				int minIndex = -1;
				double min = 0;
				double degree = pfunctionValuesMatrix.getRowVector(i).cosine(new ArrayRealVector(lambdaVectors[i]));
				double theta = k * problem.getNumberOfObjectives() * cosineLambda[i] * Math.acos(degree);
				boolean firstFlag = true;
				for (int j = 1; j < neighborhood[i].length; j++) {
					if (subRegion[neighborhood[i][j]].size() != 0) {
						if (firstFlag) {
							min = pbi(population.get((Integer) subRegion[neighborhood[i][j]].get(0)), lambdaVectors[i], theta);
							minIndex = (Integer) subRegion[neighborhood[i][j]].get(0);
							firstFlag = false;
						}
						for (int k = 0; k < subRegion[i].size(); k++) {
							double temp = pbi(population.get((Integer) subRegion[neighborhood[i][j]].get(k)), lambdaVectors[i], theta);
							if (temp < min) {
								minIndex = (int) subRegion[neighborhood[i][j]].get(k);
							}
						}
					}
					if (minIndex != -1) {
						pbestIndex[i] = minIndex;
					} else {
						pbestIndex[i] = PseudoRandom.randInt(0, archive.size() - 1);
					}
				}
			}
			//3）不止一个那么就用自适应权值来比较哪一个好
			else {
				double degree = pfunctionValuesMatrix.getRowVector(i).cosine(new ArrayRealVector(lambdaVectors[i]));
				double theta = k * problem.getNumberOfObjectives() * cosineLambda[i] * Math.acos(degree);
				double min = pbi(population.get((Integer) subRegion[i].get(0)), lambdaVectors[i], theta);
				for (int j = 1; j < subRegion[i].size(); j++) {
					double temp = pbi(population.get((Integer) subRegion[i].get(j)), lambdaVectors[i], theta);
					if (temp < min) {
						pbestIndex[i] = (int) subRegion[i].get(j);
					}
				}
			}
		}
		return pbestIndex;
	}

	private void updatePopulationPso(int[] pbestIndex) throws JMException {
		double w;
		int ran;
		double c1, c2, r1, r2;
		Variable[] pbest, gbest;
		for (int i = 0; i < populationSize; i++) {
			Variable[] paritcle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			pbest = archive.get(pbestIndex[i]).getDecisionVariables();
			ran = PseudoRandom.randInt(0, archive.size() - 1);
			gbest = archive.get(ran).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				w = PseudoRandom.randDouble(0.1, 0.5);
				r1 = PseudoRandom.randDouble();
				c1 = PseudoRandom.randDouble(1.5, 2.0);
				r2 = PseudoRandom.randDouble();
				c2 = PseudoRandom.randDouble(1.5, 2.0);
				double temp = (w * velocity[j]) + c1 * r1 * (pbest[j].getValue() - paritcle[j].getValue())
						+ c2 * r2 * (gbest[j].getValue() - paritcle[j].getValue());
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

	private void referenceSelect() {
		//Objective value Translation
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(tempPopulation.writeObjectivesToMatrix());
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
		}
		for (int i = 0; i < tempPopulation.size(); i++) {
			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
		}
		//Population Partition
		//Calculate the cosine between referenceVector and X
		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
		double[][] cosine = new double[tempPopulation.size()][populationSize];
		for (int i = 0; i < tempPopulation.size(); i++) {
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
		for (int i = 0; i < tempPopulation.size(); i++) {
			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
			subRegion[maxIndex].add(i);
		}
		//Angle-Penalized Distance(APD) Calculation
		archive.clear();
		double runIteration = Math.pow(((double) iteration / maxIterations), 2);
		double theta;
		double norm;
		for (int j = 0; j < populationSize; j++) {
			if (0 != subRegion[j].size()) {
				double[] apd = new double[subRegion[j].size()];
				for (int i = 0; i < subRegion[j].size(); i++) {
					theta = Math.acos(cosine[(int) subRegion[j].get(i)][j]) / cosineLambda[j];
					norm = functionValueMatrix.getRowVector((Integer) subRegion[j].get(i)).getNorm();
					apd[i] = (1 + problem.getNumberOfObjectives() * runIteration *
							theta) * norm;
				}
				int min = new ArrayRealVector(apd).getMinIndex();
				archive.add(tempPopulation.get((int) subRegion[j].get(min)));
			}
		}
	}

	private void offspringCreation() throws JMException {
		SolutionSet copySolution = new SolutionSet(populationSize);
		tempPopulation.clear();
		for (int j = 0; j < archive.size(); j++) {
			copySolution.add(new Solution(archive.get(j)));
			tempPopulation.add(new Solution(archive.get(j)));
		}
		for (int i = 0; i < archive.size() / 2; i++) {
			Solution[] parents = new Solution[2];
			int index = PseudoRandom.randInt(0, archive.size() - 1);
			parents[0] = copySolution.get(index);
			index = PseudoRandom.randInt(0, archive.size() - 1);
			parents[1] = copySolution.get(index);
			Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
			mutationOperator.execute(offSpring[0]);
			problem.evaluate(offSpring[0]);
			tempPopulation.add(offSpring[0]);
			problem.evaluate(offSpring[1]);
			tempPopulation.add(offSpring[1]);
		}
		if (archive.size() % 2 == 1) {
			int index;
			Solution[] parents = new Solution[2];
			index = PseudoRandom.randInt(0, archive.size() - 1);
			parents[0] = copySolution.get(index);
			index = PseudoRandom.randInt(0, archive.size() - 1);
			parents[1] = copySolution.get(index);
			Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
			mutationOperator.execute(offSpring[0]);
			problem.evaluate(offSpring[0]);
			tempPopulation.add(offSpring[0]);
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
			population.add(newSolution);
			archive.add(newSolution);
		}
	}

	private double pbi(Solution indiv, double[] lambda, double theta) {
		int i;
		double d1, d2, nl;
		double fin;
		theta = 5;

		d1 = d2 = nl = 0.0;
		for (i = 0; i < problem.getNumberOfObjectives(); i++) {
			d1 += (indiv.getObjective(i) - idealPoint[i]) * lambda[i];
			nl += Math.pow(lambda[i], 2.0);
		}
		d1 = Math.abs(d1) / Math.sqrt(nl);
		if (nl == 0.0) {
			System.out
					.println("ERROR: dived by zero(bad weihgted vector)\n");
			System.exit(0);
		}
		for (i = 0; i < problem.getNumberOfObjectives(); i++) {
			d2 += Math.pow((indiv.getObjective(i) - idealPoint[i])
					- (d1 * lambda[i]), 2.0);
		}
		d2 = Math.sqrt(d2);
		fin = (d1 + theta * d2);
		return fin;
	}

}