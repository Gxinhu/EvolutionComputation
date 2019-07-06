//package jmetal.metaheuristics.ragpso;
//
//import jmetal.core.*;
//import jmetal.qualityIndicator.QualityIndicator;
//import jmetal.util.JMException;
//import jmetal.util.NonDominatedSolutionList;
//import jmetal.util.PseudoRandom;
//import jmetal.util.createWeight;
//import jmetal.util.deepcopy.deepCopy;
//import org.apache.commons.math3.linear.Array2DRowRealMatrix;
//import org.apache.commons.math3.linear.ArrayRealVector;
//import org.apache.commons.math3.linear.RealMatrix;
//import org.apache.commons.math3.linear.RealVector;
//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Random;
//
//public class ragmopsoVersion4_2DeletePbi extends Algorithm {
//	private static final long serialVersionUID = 2107684627645440737L;
//	private Problem problem;
//	private QualityIndicator indicator;
//	private int t;
//	private int[][] neighborhood;
//	private int populationSize;
//	/**
//	 * Stores the SlutionSett
//	 */
//	private SolutionSet archive;
//	private SolutionSet population;
//	private SolutionSet tempPopulation;
//	/**
//	 * Z vector (ideal point)
//	 */
//	private double[] idealPoint;
//	private double[] nadirPoint;
//	/**
//	 * Lambda vectors
//	 */
//	private double[][] lambdaVectors;
//	private double[][] lambdaVectors0;
//	private double[] cosineLambda;
//	private Random random = new Random();
//	private int runtime;
//	private int iteration, iteration1, iteration2;
//	private Operator mutationOperator;
//	private Operator crossoverDeOperator;
//	private Operator crossoverSbxOperator;
//	private Operator selectionOperator;
//	private double w;
//	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;
//	private int maxIterations;
//
//	ragmopsoVersion4_2DeletePbi(Problem problem, QualityIndicator indicator, int i) {
//		super(problem);
//		this.problem = problem;
//		this.indicator = indicator;
//		this.runtime = i;
//	}
//
//	@Override
//	public SolutionSet execute() throws JMException, ClassNotFoundException {
//		vMax = 25.0;
//		vMin = 5.0;
//		pMax = 0.8;
//		pMin = 0.1;
//		fMax = 25.0;
//		fMin = 0.25;
//		iteration1 = maxIterations / 5;
//		iteration2 = 4 * maxIterations / 5;
//		maxIterations = (Integer) this.getInputParameter("maxIterations");
//		populationSize = (Integer) this.getInputParameter("swarmSize");
//		iteration = 0;
//		cosineLambda = new double[populationSize];
//		archive = new SolutionSet(populationSize);
//		population = new SolutionSet(populationSize);
//		tempPopulation = new SolutionSet(2 * populationSize);
//		mutationOperator = operators_.get("mutation");
//		crossoverDeOperator = operators_.get("crossoverDe");
//		crossoverSbxOperator = operators_.get("crossoverSbx");
//		selectionOperator = operators_.get("selection");
//		t = populationSize / 10;
//		neighborhood = new int[populationSize][t];
//		idealPoint = new double[problem.getNumberOfObjectives()];
//		nadirPoint = new double[problem.getNumberOfObjectives()];
//		lambdaVectors = new double[populationSize][problem
//				.getNumberOfObjectives()];
//		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
//		lambdaVectors0 = deepCopy.deepCopysDouble2d(lambdaVectors);
//		initNeighborhood();
//		initPopulation();
//		while (iteration < maxIterations) {
//			offspringCreation();
//			referenceSelect();
//			weightVectorAdaption();
//			++iteration;
//			calculateCoefficientValues();
//			offspringcreationbypso();
//			referenceselectpso();
//			weightVectorAdaption();
//			++iteration;
//		}
//		return archive;
//	}
//
//	private void calculateCoefficientValues() {
//		//calculate Vc
//		double vc;
//		if (iteration < iteration1) {
//			vc = vMax;
//		} else if (iteration > iteration1 & iteration < iteration2) {
//			vc = (iteration - iteration1) * (vMin - vMax) / (iteration2 - iteration1) + vMax;
//		} else {
//			vc = vMin;
//		}
//		//calculate p1
//		double p1;
//		if (iteration < iteration1) {
//			p1 = pMin;
//		} else if (iteration > iteration1 & iteration < (iteration2 - iteration1) / 2.0) {
//			p1 = (iteration - iteration1) * (pMax - pMin) / ((double) (iteration2 - iteration1) / 2.0 - iteration1) + pMin;
//		} else if (iteration < iteration2 & iteration > (iteration2 - iteration1) / 2.0) {
//			p1 = (iteration - iteration1) * (pMin - pMax) / (iteration2 - (double) (iteration2 - iteration1) / 2.0) + pMax;
//		} else {
//			p1 = pMin;
//		}
//		//calculate F
//		double f;
//		if (iteration < iteration1) {
//			f = fMin;
//		} else if (iteration > iteration1 & iteration < iteration2) {
//			f = 1;
//		} else {
//			f = fMax;
//		}
//		double a = Math.sqrt(f);
//		double m1 = (a + 1.0) * (a + 1.0) * (a * a + 3.0 * a + 1.0);
//		double m2 = (a + 1.0) * (a + 1.0) * (2.0 * a * a + 3.0 * a + 2.0);
//		w = (m1 * vc + m2 * p1 * vc + p1 - 1.0) / (m2 * vc + m1 * p1 * vc - p1 + 1.0);
//		double c = 2.0 * (1 - p1) * (w + 1.0) / (a + 1.0);
//		mu1 = c / 2.0;
//		mu2 = a * c / 2.0;
//		sigma1 = c / Math.sqrt(12.0);
//		sigma2 = a * c / Math.sqrt(12.0);
//	}
//
//	private void referenceselectpso() {
//		SolutionSet tempAchieve = new SolutionSet(populationSize);
//		//Objective value Translation
//		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
//		RealMatrix pfuncitionValuesMatrix = new Array2DRowRealMatrix(population.writeObjectivesToMatrix());
//		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
//			double minValue1 = functionValueMatrix.getColumnVector(i).getMinValue();
//			double minValue2 = functionValueMatrix.getColumnVector(i).getMinValue();
//			if (minValue1 < minValue2) {
//				idealPoint[i] = minValue1;
//			} else {
//				idealPoint[i] = minValue2;
//			}
//		}
//		for (int i = 0; i < archive.size(); i++) {
//			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
//		}
//		for (int i = 0; i < population.size(); i++) {
//			pfuncitionValuesMatrix.setRowVector(i, pfuncitionValuesMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
//		}
//		//Population Partition
//		//Calculate the cosine between referenceVector and X
//		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
//		double[][] cosine = new double[archive.size()][populationSize];
//		for (int i = 0; i < archive.size(); i++) {
//			for (int j = 0; j < populationSize; j++) {
//				RealVector functionValueVector = functionValueMatrix.getRowVector(i);
//				RealVector referenceVector = referenceMatrix.getRowVector(j);
//				cosine[i][j] = functionValueVector.cosine(referenceVector);
//			}
//		}
//		//Partition
//		List[] subRegion = new List[populationSize];
//		for (int i = 0; i < populationSize; i++) {
//			subRegion[i] = new ArrayList();
//		}
//		RealMatrix cosineMatrix = new Array2DRowRealMatrix(cosine);
//		int maxIndex;
//		for (int i = 0; i < archive.size(); i++) {
//			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
//			subRegion[maxIndex].add(i);
//		}
//		//Angle-Penalized Distance(APD) Calculation
//		for (int i = 0; i < populationSize; i++) {
//			if (subRegion[i].size() != 0) {
//				double degree = pfuncitionValuesMatrix.getRowVector(i).cosine(new ArrayRealVector(lambdaVectors[i]));
//				double theta = k * problem.getNumberOfObjectives() * cosineLambda[i] * Math.acos(degree);
//				tempAchieve.add(chooseSolution(population.get(i), subRegion[i], theta, i));
//			} else {
//				tempAchieve.add(new Solution(population.get(i)));
//			}
//		}
//		archive = tempAchieve;
//	}
//
//	private Solution chooseSolution(Solution solution, List list, double theta, int index) {
//		NonDominatedSolutionList tempSolutions = new NonDominatedSolutionList(new jmetal.util.comparators.DominanceComparator());
//		tempSolutions.add(solution);
//		for (Object o : list) {
//			tempSolutions.add(archive.get((Integer) o));
//		}
//		if (tempSolutions.size() == 1) {
//			return tempSolutions.get(0);
//		} else {
//			double minpbi1 = pbi(tempSolutions.get(0), lambdaVectors[index], theta);
//			int minindex = 0;
//			for (int j = 1; j < tempSolutions.size(); j++) {
//				double minpbi2 = pbi(tempSolutions.get(j), lambdaVectors[index], theta);
//				if (minpbi2 < minpbi1) {
//					minpbi1 = minpbi2;
//					minindex = j;
//				}
//			}
//			return tempSolutions.get(minindex);
//		}
//	}
//
//	private void offspringcreationbypso() throws JMException {
//		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
//		RealMatrix pfunctionValuesMatrix = new Array2DRowRealMatrix(population.writeObjectivesToMatrix());
//		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
//			idealPoint[i] = (functionValueMatrix.getColumnVector(i).getMinValue() < pfunctionValuesMatrix.getColumnVector(i).getMinValue())
//					? functionValueMatrix.getColumnVector(i).getMinValue() : pfunctionValuesMatrix.getColumnVector(i).getMinValue();
//
//		}
//		for (int i = 0; i < archive.size(); i++) {
//			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
//			pfunctionValuesMatrix.setRowVector(i, pfunctionValuesMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
//		}
//		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
//		double[][] cosine = new double[archive.size()][populationSize];
//		for (int i = 0; i < archive.size(); i++) {
//			for (int j = 0; j < populationSize; j++) {
//				RealVector functionValueVector = functionValueMatrix.getRowVector(i);
//				RealVector referenceVector = referenceMatrix.getRowVector(j);
//				cosine[i][j] = functionValueVector.cosine(referenceVector);
//			}
//		}
//		//Partition
//		List[] subRegion = new List[populationSize];
//		for (int i = 0; i < populationSize; i++) {
//			subRegion[i] = new ArrayList();
//		}
//		RealMatrix cosineMatrix = new Array2DRowRealMatrix(cosine);
//		int maxIndex;
//		double[] angle = new double[populationSize];
//		for (int i = 0; i < archive.size(); i++) {
//			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
//			angle[i] = cosineMatrix.getRowVector(i).getMaxValue();
//			subRegion[maxIndex].add(i);
//		}
//		int[] pbestIndex = findPbest(subRegion, functionValueMatrix);
//		updatePopulationPso(pbestIndex);
//
//	}
//
//	private int[] findPbest(List[] subRegion, RealMatrix functionValuesMatrix) {
//		int[] pbestIndex = new int[populationSize];
//		for (int i = 0; i < populationSize; i++) {
//			//1)当前权值向量为空那个么在周围邻居里面寻找一个作为pbset 当然要是周围寻找的还是为空的话就通过邻居的远近距离找个最近的，或者在所有邻居里面找然后通过自适应pbi来寻找pbest
//			if (subRegion[i].size() == 1) {
//				pbestIndex[i] = (int) subRegion[i].get(0);
//			}
//			//2）当前去去权值向量里面的subregion里面就只有一个个体那么就直接选择这个个体作为pbest
//			else if (subRegion[i].size() == 0) {
//				int minIndex = -1;
//				double min = 0;
//				double degree = functionValuesMatrix.getRowVector(i).cosine(new ArrayRealVector(lambdaVectors[i]));
//				double theta = Math.acos(degree / cosineLambda[j]);
//				double norm = functionValuesMatrix.getRowVector((Integer) subRegion[j].get(i)).getNorm();
//				double runIteration=Math.pow(iteration/maxIterations,2);
//				apd[i] = (1 + problem.getNumberOfObjectives() * runIteration *
//						theta) * norm;
//				boolean firstFlag = true;
//				for (int j = 1; j < neighborhood[i].length; j++) {
//					if (subRegion[neighborhood[i][j]].size() != 0) {
//						if (firstFlag) {
//							min = pbi(population.get((Integer) subRegion[neighborhood[i][j]].get(0)), lambdaVectors[i], theta);
//							minIndex = (Integer) subRegion[neighborhood[i][j]].get(0);
//							firstFlag = false;
//						}
//						for (int k = 0; k < subRegion[i].size(); k++) {
//							double temp = pbi(population.get((Integer) subRegion[neighborhood[i][j]].get(k)), lambdaVectors[i], theta);
//							if (temp < min) {
//								minIndex = (int) subRegion[neighborhood[i][j]].get(k);
//							}
//						}
//					}
//					if (minIndex != -1) {
//						pbestIndex[i] = minIndex;
//					} else {
//						pbestIndex[i] = PseudoRandom.randInt(0, archive.size() - 1);
//					}
//				}
//			}
//			//3）不止一个那么就用自适应权值来比较哪一个好
//			else {
//				double degree = functionValuesMatrix.getRowVector((Integer) subRegion[i].get(0)).cosine(new ArrayRealVector(lambdaVectors[i]));
//				double theta = Math.acos(degree / cosineLambda[j]);
//				double norm = functionValuesMatrix.getRowVector((Integer) subRegion[j].get(i)).getNorm();
//				double runIteration=Math.pow(iteration/maxIterations,2);
//				apd[i] = (1 + problem.getNumberOfObjectives() * runIteration *
//						theta) * norm;
//				double min = pbi(population.get((Integer) subRegion[i].get(0)), lambdaVectors[i], theta);
//				for (int j = 1; j < subRegion[i].size(); j++) {
//					double temp = pbi(population.get((Integer) subRegion[i].get(j)), lambdaVectors[i], theta);
//					if (temp < min) {
//						pbestIndex[i] = (int) subRegion[i].get(j);
//					}
//				}
//			}
//		}
//		return pbestIndex;
//	}
//
//	private void updatePopulationPso(int[] pbestIndex) throws JMException {
//		int ran;
//		Variable[] pbest, gbest;
//		for (int i = 0; i < populationSize; i++) {
//			Variable[] paritcle = population.get(i).getDecisionVariables();
//			double[] velocity = population.get(i).getSpeed();
//			pbest = archive.get(pbestIndex[i]).getDecisionVariables();
//			ran = PseudoRandom.randInt(0, archive.size() - 1);
//			gbest = archive.get(ran).getDecisionVariables();
//			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
//				double theta1 = random.nextGaussian() * sigma1 + mu1;
//				double theta2 = random.nextGaussian() * sigma2 + mu2;
//				double temp = (w * velocity[j]) + theta1 * (pbest[j].getValue() - paritcle[j].getValue())
//						+ theta2 * (gbest[j].getValue() - paritcle[j].getValue());
//				population.get(i).setSpeed(j, temp);
//			}
//		}
//		for (int n = 0; n < this.populationSize; n++) {
//
//			// DecisionVariables particle = ;
//			for (int var = 0; var < population.get(n).getDecisionVariables().length; var++) {
//				population.get(n).getDecisionVariables()[var].setValue(population.get(n).getDecisionVariables()[var].getValue() + population.get(n).getSpeed()[var]);
//				if (population.get(n).getDecisionVariables()[var].getValue() < problem.getLowerLimit(var)) {
//					population.get(n).getDecisionVariables()[var].setValue(problem.getLowerLimit(var));
//					population.get(n).setSpeed(var, 0.0);
//				}
//				if (population.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
//					population.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
//					population.get(n).setSpeed(var, 0.0);
//				}
//			}
//			problem.evaluate(population.get(n));
//		}
//	}
//
//	private void referenceSelect() {
//		//Objective value Translation
//		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(tempPopulation.writeObjectivesToMatrix());
//		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
//			idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
//		}
//		for (int i = 0; i < tempPopulation.size(); i++) {
//			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
//		}
//		//Population Partition
//		//Calculate the cosine between referenceVector and X
//		RealMatrix referenceMatrix = new Array2DRowRealMatrix(lambdaVectors);
//		double[][] cosine = new double[tempPopulation.size()][populationSize];
//		for (int i = 0; i < tempPopulation.size(); i++) {
//			for (int j = 0; j < populationSize; j++) {
//				RealVector functionValueVector = functionValueMatrix.getRowVector(i);
//				RealVector referenceVector = referenceMatrix.getRowVector(j);
//				cosine[i][j] = functionValueVector.cosine(referenceVector);
//			}
//		}
//		//Partition
//		List[] subRegion = new List[populationSize];
//		for (int i = 0; i < populationSize; i++) {
//			subRegion[i] = new ArrayList();
//		}
//		RealMatrix cosineMatrix = new Array2DRowRealMatrix(cosine);
//		int maxIndex;
//		for (int i = 0; i < tempPopulation.size(); i++) {
//			maxIndex = cosineMatrix.getRowVector(i).getMaxIndex();
//			subRegion[maxIndex].add(i);
//		}
//		//Angle-Penalized Distance(APD) Calculation
//		archive.clear();
//		double runIteration = Math.pow(((double) iteration / maxIterations), 2);
//		double theta;
//		double norm;
//		for (int j = 0; j < populationSize; j++) {
//			if (0 != subRegion[j].size()) {
//				double[] apd = new double[subRegion[j].size()];
//				for (int i = 0; i < subRegion[j].size(); i++) {
//					theta = Math.acos(cosine[(int) subRegion[j].get(i)][j]) / cosineLambda[j];
//					norm = functionValueMatrix.getRowVector((Integer) subRegion[j].get(i)).getNorm();
//					apd[i] = (1 + problem.getNumberOfObjectives() * runIteration *
//							theta) * norm;
//				}
//				int min = new ArrayRealVector(apd).getMinIndex();
//				archive.add(tempPopulation.get((int) subRegion[j].get(min)));
//			}
//		}
//	}
//
//	private void offspringCreation() throws JMException {
//		SolutionSet copySolution = new SolutionSet(populationSize);
//		tempPopulation.clear();
//		for (int j = 0; j < archive.size(); j++) {
//			copySolution.add(new Solution(archive.get(j)));
//			tempPopulation.add(new Solution(archive.get(j)));
//		}
//		for (int i = 0; i < archive.size(); i++) {
//			int k = PseudoRandom.randInt(1, 3);
//			Solution offSpring;
//			if (k == 1) {
//				Solution parents[] = (Solution[]) selectionOperator.execute(new Object[]{
//						archive, i});
//				offSpring = (Solution) crossoverDeOperator.execute(new Object[]{archive.get(i), parents});
//			} else if (k == 2) {
//				int s = PseudoRandom.randInt(0, archive.size());
//				Solution[] praents = new Solution[2];
//				praents[0] = archive.get(i);
//				praents[1] = archive.get(PseudoRandom.randInt(0, archive.size() - 1));
//				Solution[] offSprings = (Solution[]) crossoverSbxOperator.execute(praents);
//				offSpring = offSprings[0];
//			} else {
//				offSpring = (Solution) mutationOperator.execute(archive.get(i));
//			}
//
//			problem.evaluate(offSpring);
//			tempPopulation.add(offSpring);
//		}
//	}
//
//	private void weightVectorAdaption() {
//		if (iteration % Math.ceil(maxIterations * 0.1) == 0) {
//			RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
//			RealVector subtract = new ArrayRealVector(problem.getNumberOfObjectives());
//			for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
//				idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
//				nadirPoint[i] = functionValueMatrix.getColumnVector(i).getMaxValue();
//				subtract.setEntry(i, nadirPoint[i] - idealPoint[i]);
//			}
//			RealMatrix lambdaMatrix = new Array2DRowRealMatrix(lambdaVectors0);
//			for (int i = 0; i < populationSize; i++) {
//				lambdaMatrix.setRowVector(i, lambdaMatrix.getRowVector(i).ebeMultiply(subtract));
//				lambdaMatrix.setRowVector(i, lambdaMatrix.getRowVector(i).mapDivide(lambdaMatrix.getRowVector(i).getNorm()));
//			}
//			this.lambdaVectors = lambdaMatrix.getData();
//			this.initNeighborhood();
//		}
//
//
//	}
//
//	public void initNeighborhood() {
//		RealMatrix lambdaMatrix = new Array2DRowRealMatrix(lambdaVectors);
//		lambdaMatrix = lambdaMatrix.multiply(lambdaMatrix.transpose());
//		for (int k = 0; k < populationSize; k++) {
//			double[] arrays = lambdaMatrix.getRowVector(k).toArray();
//			ArrayList<Integer> index = new ArrayList<>(arrays.length);
//			for (int i = 0; i < arrays.length; i++) {
//				index.add(i);
//			}
//			index.sort((o1, o2) -> Double.compare(arrays[o2], arrays[o1]));
//			for (int j = 0; j < populationSize; j++) {
//				if (j < t) {
//					neighborhood[k][j] = index.get(j);
//				}
//			}
//			cosineLambda[k] = Math.acos(arrays[index.get(1)]);
//		}
//	}
//
//	public void initPopulation() throws JMException,
//			ClassNotFoundException {
//		for (int i = 0; i < populationSize; i++) {
//			Solution newSolution = new Solution(problem);
//			problem.evaluate(newSolution);
//			if (this.problem.getNumberOfConstraints() != 0) {
//				problem.evaluateConstraints(newSolution);
//			}
//			// evaluations++;
//			population.add(new Solution(newSolution));
//			archive.add(new Solution(newSolution));
//		}
//	}
//
//
//}