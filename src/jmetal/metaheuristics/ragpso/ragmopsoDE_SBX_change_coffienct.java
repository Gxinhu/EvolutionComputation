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

public class ragmopsoDE_SBX_change_coffienct extends Algorithm {
	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private QualityIndicator indicator;
	private int t;
	private int[][] neighborhood;
	private int populationSize;
	private double judge;
	private double equal;
	/**
	 * Stores the SolutionSett
	 */
	private SolutionSet archive;
	private SolutionSet population;
	private SolutionSet tempPopulation;
	private SolutionSet lastPbest;
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
	private Operator mutationOperator;
	private Operator crossoverDeOperator;
	private Operator crossoverSbxOperator;
	private Operator selectionOperator;
	private double w;
	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;
	private int maxIterations;

	ragmopsoDE_SBX_change_coffienct(Problem problem, QualityIndicator indicator, int i) {
		super(problem);
		this.problem = problem;
		this.indicator = indicator;
		this.runtime = i;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		judge = 0;
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
		mutationOperator = operators_.get("mutation");
		crossoverDeOperator = operators_.get("crossoverDe");
		crossoverSbxOperator = operators_.get("crossoverSbx");
		selectionOperator = operators_.get("selection");
		iteration = 0;
		cosineLambda = new double[populationSize];
		archive = new SolutionSet(populationSize);
		population = new SolutionSet(populationSize);
		lastPbest = new SolutionSet(populationSize);
		tempPopulation = new SolutionSet(2 * populationSize);
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
			offspringCreationbypso();
			referenceselectpso();
			weightVectorAdaption();
			++iteration;
		}
		return population;
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
		double condition = judge / populationSize;
		for (int i = 0; i < archive.size(); i++) {
			double k = PseudoRandom.randDouble();
			Solution offSpring;
			if (k < condition) {
				Solution[] parents = (Solution[]) selectionOperator.execute(new Object[]{
						archive, i});
				offSpring = (Solution) crossoverDeOperator.execute(new Object[]{archive.get(i), parents});
			} else {
				Solution[] parents = new Solution[2];
				parents[0] = archive.get(i);
				parents[1] = archive.get(PseudoRandom.randInt(0, archive.size() - 1));
				Solution[] offSprings = (Solution[]) crossoverSbxOperator.execute(parents);
				offSpring = offSprings[0];
			}
			mutationOperator.execute(offSpring);


			problem.evaluate(offSpring);
			tempPopulation.add(offSpring);
		}
	}

	private void calculateCoefficientValues(int index, RealVector solution) {
		//calculate Vc
		double distance = 0;
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			distance += nadirPoint[i] - idealPoint[i];
		}
		double vc = distance / problem.getNumberOfObjectives();
		//calculate p1
		double p1;
		RealVector referenceVector = new ArrayRealVector(lambdaVectors[index]);
		p1 = 1 - Math.acos(solution.cosine(referenceVector)) / cosineLambda[index];
//		if (iteration < iteration1) {
//			p1 = pMin;
//		} else if (iteration > iteration1 & iteration < (iteration2 - iteration1) / 2.0) {
//			p1 = (iteration - iteration1) * (pMax - pMin) / ((double) (iteration2 - iteration1) / 2.0 - iteration1) + pMin;
//		} else if (iteration < iteration2 & iteration > (iteration2 - iteration1) / 2.0) {
//			p1 = (iteration - iteration1) * (pMin - pMax) / (iteration2 - (double) (iteration2 - iteration1) / 2.0) + pMax;
//		} else {
//			p1 = pMin;
//		}
		//calculate F
		double f;
//		if (iteration < iteration1) {
//			f = fMin;
//		} else if (iteration > iteration1 & iteration < iteration2) {
//			f = 1;
//		} else {
//			f = fMax;
//		}
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
			nadirPoint[i] = (functionValueMatrix.getColumnVector(i).getMaxValue() < pfunctionValuesMatrix.getColumnVector(i).getMaxValue())
					? functionValueMatrix.getColumnVector(i).getMaxValue() : pfunctionValuesMatrix.getColumnVector(i).getMaxValue();

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
		int[] pBestIndex = findpbest(subRegion, cosine, functionValueMatrix, referenceMatrix);
		updatePopulationPso(pBestIndex, functionValueMatrix, pfunctionValuesMatrix);

	}

	private int[] findpbest(List[] subRegion, double[][] cosine, RealMatrix valueMatrix, RealMatrix functionValueMatrix) {
		int[] pbestIndex = new int[populationSize];
		for (int i = 0; i < populationSize; i++) {
			//1)当前权值向量为空那个么在周围邻居里面寻找一个作为pbset 当然要是周围寻找的还是为空的话就通过邻居的远近距离找个最近的，或者在所有邻居里面找然后通过自适应pbi来寻找pbest
			if (subRegion[i].size() == 1) {
				pbestIndex[i] = (int) subRegion[i].get(0);
			}
			//2）当前去去权值向量里面的subregion里面就只有一个个体那么就直接选择这个个体作为pbest
			else if (subRegion[i].size() == 0) {
				double runIteration = Math.pow((double) iteration / maxIterations, 2);
				double theta = Math.acos(cosine[0][i]) / cosineLambda[i];
				double norm = functionValueMatrix.getRowVector(0).getNorm();
				double minapd1 = (1 + problem.getNumberOfObjectives() * runIteration *
						theta) * norm;
				int minindex = 0;
				for (int j = 1; j < archive.size(); j++) {
					theta = Math.acos(cosine[j][i]) / cosineLambda[i];
					norm = functionValueMatrix.getRowVector(j).getNorm();
					double minapd2 = (1 + problem.getNumberOfObjectives() * runIteration *
							theta) * norm;
					if (minapd2 < minapd1) {
						minindex = j;
						minapd1 = minapd2;
					}
				}
				pbestIndex[i] = minindex;
			}
			//3）不止一个那么就用自适应权值来比较哪一个好
			else {
				pbestIndex[i] = chooseSolutionIndex(subRegion[i], i, cosine, functionValueMatrix);
			}
		}
//		for (int i = 0; i < populationSize; i++) {
//			if (pbestIndex[i] == -1) {
//				for (int j = 0; j < neighborhood[i].length; j++) {
//					if (pbestIndex[neighborhood[i][j]] != -1) {
//						pbestIndex[i] = pbestIndex[neighborhood[i][j]];
//						break;
//					}
//				}
//				while (pbestIndex[i] == -1) {
//					pbestIndex[i] = PseudoRandom.randInt(0, archive.size() - 1);
//				}
//			}
//		}
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

	private void updatePopulationPso(int[] pBestIndex, RealMatrix functionValueMatrix, RealMatrix pfunctionValueMatrix) throws JMException {
		int ran;
		Variable[] pbest, gbest;
		RealMatrix pbsetMatrix = new Array2DRowRealMatrix(lastPbest.writeObjectivesToMatrix());
		judge = 0;
		equal = 0;
		double c1, c2, r1, r2;
		for (int i = 0; i < populationSize; i++) {
			Variable[] particle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			//choose pbest
			RealVector referenceVector = new ArrayRealVector(lambdaVectors[i]);
			RealVector functionValueVector = functionValueMatrix.getRowVector(pBestIndex[i]);
			double apd1;
			double apd2;
			double cosine = functionValueVector.cosine(referenceVector);
			double runIteration = Math.pow((double) iteration / maxIterations, 2);
			double theta = Math.acos(cosine) / cosineLambda[i];
			double norm = functionValueVector.getNorm();
			apd1 = (1 + problem.getNumberOfObjectives() * runIteration *
					theta) * norm;
			functionValueVector = pbsetMatrix.getRowVector(i);
			cosine = functionValueVector.cosine(referenceVector);
			theta = Math.acos(cosine) / cosineLambda[i];
			norm = functionValueVector.getNorm();
			apd2 = (1 + problem.getNumberOfObjectives() * runIteration *
					theta) * norm;
			if (apd1 < apd2) {
				lastPbest.replace(i, new Solution(archive.get(pBestIndex[i])));
			} else {
				judge += 1;
			}
			pbest = lastPbest.get(i).getDecisionVariables();

			// Update by PSO

			calculateCoefficientValues(i, pfunctionValueMatrix.getRowVector(i));
			ran = PseudoRandom.randInt(0, archive.size() - 1);
			gbest = archive.get(ran).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				double theta1 = random.nextGaussian() * sigma1 + mu1;
				double theta2 = random.nextGaussian() * sigma2 + mu2;
				double temp = (w * velocity[j]) + theta1 * (pbest[j].getValue() - particle[j].getValue())
						+ theta2 * (gbest[j].getValue() - particle[j].getValue());
				population.get(i).setSpeed(j, temp);

			}
		}
		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle  ;
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
			lastPbest.add(new Solution(newSolution));
		}
	}

}