package jmetal.metaheuristics.rvea;


import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList2;
import jmetal.util.PseudoRandom;
import jmetal.util.createWeight;
import jmetal.util.deepcopy.deepCopy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class rvea extends Algorithm {

	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private QualityIndicator indicator;
	private int t;
	private int[][] neighborhood;
	private int populationSize;
	/**
	 * Stores the population
	 */
	private SolutionSet population;
	private SolutionSet tempPopulation;
	private SolutionSet nadirAchieve;
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

	public rvea(Problem problem, QualityIndicator indicator, int i) {
		super(problem);
		this.problem = problem;
		this.indicator = indicator;
		this.runtime = i;
	}

	public rvea(Problem problem, int i) {
		super(problem);
		this.problem = problem;
		this.runtime = i;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		iteration = 0;
		cosineLambda = new double[populationSize];
		population = new SolutionSet(populationSize);
		tempPopulation = new SolutionSet(2 * populationSize);
		nadirAchieve = new SolutionSet(populationSize);
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
//		estimatedPoint();
		while (iteration < maxIterations) {
			offspringCreation();
			referenceSelect();
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
		population.clear();
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
				population.add(tempPopulation.get((int) subRegion[j].get(min)));
			}
		}
	}

	private void offspringCreation() throws JMException {
		SolutionSet copySolution = new SolutionSet(populationSize);
		tempPopulation.clear();
		for (int j = 0; j < population.size(); j++) {
			copySolution.add(new Solution(population.get(j)));
			tempPopulation.add(new Solution(population.get(j)));
		}
		for (int i = 0; i < population.size() / 2; i++) {
			Solution[] parents = new Solution[2];
			int index = PseudoRandom.randInt(0, population.size() - 1);
			parents[0] = copySolution.get(index);
			index = PseudoRandom.randInt(0, population.size() - 1);
			parents[1] = copySolution.get(index);
			Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
			mutationOperator.execute(offSpring[0]);
			problem.evaluate(offSpring[0]);
			tempPopulation.add(offSpring[0]);
			problem.evaluate(offSpring[1]);
			tempPopulation.add(offSpring[1]);
		}
		if (population.size() % 2 == 1) {
			int index;
			Solution[] parents = new Solution[2];
			index = PseudoRandom.randInt(0, population.size() - 1);
			parents[0] = copySolution.get(index);
			index = PseudoRandom.randInt(0, population.size() - 1);
			parents[1] = copySolution.get(index);
			Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
			mutationOperator.execute(offSpring[0]);
			problem.evaluate(offSpring[0]);
			tempPopulation.add(offSpring[0]);
		}
	}

	private void weightVectorAdaption() {
		if (iteration % Math.ceil(maxIterations * 0.1) == 0) {
			RealMatrix functionValueMatrix = new Array2DRowRealMatrix(population.writeObjectivesToMatrix());
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
//			nadirAchieve.clear();
		}


	}

	//另一种计算
	private void estimatedPoint() {
		SolutionSet temp = new SolutionSet(populationSize + nadirAchieve.size());
		for (int i = 0; i < populationSize; i++) {
			temp.add(population.get(i));
		}
		for (int i = 0; i < this.nadirAchieve.size(); i++) {
			temp.add(this.nadirAchieve.get(i));
		}
		NonDominatedSolutionList2 list = new NonDominatedSolutionList2();
		for (int j = 0; j < temp.size(); j++) {
			list.add(temp.get(j));
		}
		RealMatrix matrix = new Array2DRowRealMatrix(list.writeObjectivesToMatrix());
		nadirAchieve.clear();
		RealVector vector;
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			vector = matrix.getColumnVector(i);
			nadirPoint[i] = vector.getMaxValue();
			nadirAchieve.add(list.get(vector.getMaxIndex()));
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
			Collections.sort(index, new Comparator<Integer>() {
				@Override
				public int compare(Integer o1, Integer o2) {
					return Double.compare(arrays[o2], arrays[o1]);
				}
			});
			for (int j = 0; j < populationSize; j++) {
				if (j < t) {
					neighborhood[k][j] = index.get(j);
				}
			}
			cosineLambda[k] = Math.acos(arrays[index.get(1)]);
		}
	} // initNeighborhood

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
		}
	} // initPopulation

} // MOPSOD