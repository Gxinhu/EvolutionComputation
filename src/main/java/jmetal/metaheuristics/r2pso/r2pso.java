package jmetal.metaheuristics.r2pso;

import jmetal.core.*;
import jmetal.metaheuristics.r2pso.util.r2calculate;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;
import jmetal.util.comparators.DegreeComparator;
import jmetal.util.deepcopy.deepCopy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class r2pso extends Algorithm {
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
	private SolutionSet archive;
	private SolutionSet population;
	private SolutionSet tempPopulation;
	private SolutionSet clonePopulation;
	private SolutionSet lastPbestPopulation;
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
	private int iteration, iteration1, iteration2;
	private Operator mutationOperator;
	private Operator crossoverOperator;
	private Operator cloneOperator;
	private double w;
	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;
	private int maxIterations;
	r2calculate R2 = new r2calculate();
	private QualityIndicator indicators;
	private Degree degree;

	r2pso(Problem problem, QualityIndicator indicator) {
		super(problem);
		this.problem = problem;
		this.indicators = indicator;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		vMax = 25.0;
		vMin = 5.0;
		pMax = 0.8;
		pMin = 0.1;
		fMax = 25.0;
		fMin = 0.25;
		degree = new Degree();
		iteration1 = maxIterations / 5;
		iteration2 = 4 * maxIterations / 5;
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
		idealPoint = new double[problem.getNumberOfObjectives()];
		nadirPoint = new double[problem.getNumberOfObjectives()];
		lambdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];
		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
		lambdaVectors0 = deepCopy.deepCopysDouble2d(lambdaVectors);
		initNeighborhood();
		initPopulation();
		while (iteration < maxIterations) {
			cloneOffspringCreation();
//			selection();
//			weightVectorAdaption();
			++iteration;
//			calculateCoefficientValues();
//			offspringcreationbypso();
//			referenceselectpso();
//			weightVectorAdaption();
//			++iteration;
		}
		return archive;
	}

	private void selection() {

	}

	private void cloneOffspringCreation() throws JMException {
		// clone operation
		int clonesize = (int) populationSize / 5;
		RealMatrix functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
		int[] extrmePoints = new int[problem.getNumberOfObjectives()];
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
			extrmePoints[i] = functionValueMatrix.getColumnVector(i).getMinIndex();
			nadirPoint[i] = functionValueMatrix.getColumnVector(i).getMaxValue() - idealPoint[i];
		}
		for (int i = 0; i < archive.size(); i++) {
			functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
		}
		//Calculate the cosine between Objective
		double cosine;
		if (archive.size() != 1) {
			for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
				RealVector functionValueVector1 = functionValueMatrix.getRowVector(extrmePoints[i]);
				for (int j = 0; j < archive.size(); j++) {
					RealVector functionValueVector2 = functionValueMatrix.getRowVector(j);
					cosine = functionValueVector1.cosine(functionValueVector2);
					if (Math.abs(cosine - 1) < 1e-8) {
						cosine = 1;
					}
					archive.get(j).setCosine(Math.acos(cosine), i);
				}
			}
		}
		degree.crowdingDegreeAssignment(archive, problem.getNumberOfObjectives());
		archive.sort(new DegreeComparator());
		SolutionSet cpopulation = new SolutionSet(clonesize);
		for (int k = 0; k < archive.size() && k < clonesize; k++) {
			cpopulation.add(new Solution(archive.get(k)));
		} // for
		clonePopulation = (SolutionSet) cloneOperator.execute(cpopulation);
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
		for (int k = 0; k < tempPopulation.size(); k++) {
			archive.add(tempPopulation.get(k));
			if (archive.size() == populationSize + 1) {
				functionValueMatrix = new Array2DRowRealMatrix(archive.writeObjectivesToMatrix());
				extrmePoints = new int[problem.getNumberOfObjectives()];
				for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
					idealPoint[i] = functionValueMatrix.getColumnVector(i).getMinValue();
					extrmePoints[i] = functionValueMatrix.getColumnVector(i).getMinIndex();
					nadirPoint[i] = functionValueMatrix.getColumnVector(i).getMaxValue() - idealPoint[i];
				}
				for (int i = 0; i < archive.size(); i++) {
					functionValueMatrix.setRowVector(i, functionValueMatrix.getRowVector(i).subtract(new ArrayRealVector(idealPoint)));
				}
				//Calculate the cosine between Objective

				for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
					RealVector functionValueVector1 = functionValueMatrix.getRowVector(extrmePoints[i]);
					for (int j = 0; j < archive.size(); j++) {
						RealVector functionValueVector2 = functionValueMatrix.getRowVector(j);
						cosine = functionValueVector1.cosine(functionValueVector2);
						if (Math.abs(cosine - 1) < 1e-8) {
							cosine = 1;
						}
						archive.get(j).setCosine(Math.acos(cosine), i);
					}
				}
				degree.crowdingDegreeAssignment(archive, problem.getNumberOfObjectives());
				archive.sort(new DegreeComparator());
				archive.remove(populationSize);
			}
		}
//		for (int i = 0; i < archive.size(); i++) {
//			for (int j = 0; j < archive.size(); j++) {
//				if (i != j) {
//					RealVector functionValueVector1 = functionValueMatrix.getRowVector(i);
//					RealVector functionValueVector2 = functionValueMatrix.getRowVector(j);
//					cosine[i][j] = functionValueVector1.cosine(functionValueVector2);
//				}
//			    else{
//			    	cosine[i][j]=-1.0;
//				}
//			}
//		}
//		for (int j = 0; j < archive.size(); j++) {
//				RealVector nadValueVector = new ArrayRealVector(nadirPoint);
//				RealVector functionValueVector2 = functionValueMatrix.getRowVector(j);
//				cosine[archive.size()][j] = nadValueVector.cosine(functionValueVector2);
//
//		}
//		RealMatrix cosineMatrix=new Array2DRowRealMatrix(cosine);
//		for(int i=0;i<=archive.size();i++){
//			if(i==archive.size()){
//				int index=cosineMatrix.getRowVector(i).getMaxIndex();
//				if(cosineMatrix.getEntry(archive.size(),index)<archive.get(index).getCosine()){
//					archive.get(index).setCosine(cosineMatrix.getEntry(archive.size(),index));
//				}
//				break;
//			}
//			archive.get(i).setCloseindex(cosineMatrix.getRowVector(i).getMaxIndex());
//			archive.get(i).setCosine(Math.acos(cosineMatrix.getRowVector(i).getEntry(archive.get(i).getCloseindex())));
//		}
//		archive.sort(new DegreeComparator());
//		clonePopulation.clear();
//		double sumDegree=0;
//		for(int i=0;i<archive.size();i++){
//			sumDegree+=archive.get(i).getCosine();
//		}
//		for(int i=0;i<archive.size()&i<t;i++){
//			for(int j=0;j<Math.ceil(populationSize*archive.get(i).getCosine()/sumDegree);j++){
//				clonePopulation.add(new Solution(archive.get(i)));
//			}
//		}
//		// operation;
//		for (int i = 0; i < clonePopulation.size(); i++) {
//			Solution[] particle2 = new Solution[2];
//			int ran;
//			particle2[0] = clonePopulation.get(i);
//			ran = PseudoRandom.randInt(0, clonePopulation.size() - 1);
//			particle2[1] = clonePopulation.get(ran);
//			Solution[] offSpring = (Solution[]) crossoverOperator.execute(particle2);
//			mutationOperator.execute(offSpring[0]);
//			problem.evaluate(offSpring[0]);
//			archive.add(offSpring[0]);
//		}

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
			p1 = pMin;
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
				tempAchieve.add(chooseSolution(population.get(i), subRegion[i], theta, i));
			} else {
				tempAchieve.add(new Solution(population.get(i)));
			}
		}
		archive = tempAchieve;
	}

	private Solution chooseSolution(Solution solution, List list, double theta, int index) {
		NonDominatedSolutionList tempSolutions = new NonDominatedSolutionList(new jmetal.util.comparators.DominanceComparator());
		tempSolutions.add(solution);
		for (Object o : list) {
			tempSolutions.add(archive.get((Integer) o));
		}
		if (tempSolutions.size() == 1) {
			return tempSolutions.get(0);
		} else {
			double minpbi1 = pbi(tempSolutions.get(0), lambdaVectors[index], theta);
			int minindex = 0;
			for (int j = 1; j < tempSolutions.size(); j++) {
				double minpbi2 = pbi(tempSolutions.get(j), lambdaVectors[index], theta);
				if (minpbi2 < minpbi1) {
					minpbi1 = minpbi2;
					minindex = j;
				}
			}
			return tempSolutions.get(minindex);
		}
	}

	private void offspringcreationbypso() throws JMException {
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
		int ran;
		Variable[] pbest, gbest;
		for (int i = 0; i < populationSize; i++) {
			Variable[] paritcle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			pbest = archive.get(pbestIndex[i]).getDecisionVariables();
			ran = PseudoRandom.randInt(0, archive.size() - 1);
			gbest = archive.get(ran).getDecisionVariables();
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
	}

	private void offspringCreation() throws JMException {
		SolutionSet copySolution = new SolutionSet(populationSize);
		tempPopulation.clear();
		for (int j = 0; j < archive.size(); j++) {
			copySolution.add(new Solution(archive.get(j)));
			tempPopulation.add(new Solution(archive.get(j)));
		}
		for (int i = 0; i < archive.size(); i++) {
			int k = PseudoRandom.randInt(1, 3);
			Solution offSpring;
			int s = PseudoRandom.randInt(0, archive.size());
			Solution[] parents = new Solution[2];
			parents[0] = archive.get(i);
			parents[1] = archive.get(PseudoRandom.randInt(0, archive.size() - 1));
			Solution[] offSprings = (Solution[]) crossoverOperator.execute(parents);
			offSpring = offSprings[0];
			offSpring = (Solution) mutationOperator.execute(archive.get(i));


			problem.evaluate(offSpring);
			tempPopulation.add(offSpring);
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