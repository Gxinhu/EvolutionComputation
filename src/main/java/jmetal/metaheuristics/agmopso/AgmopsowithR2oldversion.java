package jmetal.metaheuristics.agmopso;

import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.createWeight;
import jmetal.util.deepcopy.deepCopy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

public class AgmopsowithR2oldversion extends Algorithm {
	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;

	int run;
	int T_;
	int[][] neighborhood_;
	private double[][] lambdaVectors0;
	private double[] nadirPoint;
	public String curDir = System.getProperty("user.dir");
	/**
	 * Stores the population

	 */
	private int populationSize;
	private SolutionSet population, temppopulation, cpopulation, leader_ind;
	/**
	 * Z vector (ideal point)
	 */
	double[] idealPoint;
	/**
	 * Lambda vectors
	 */
	private double max_d = Double.MIN_VALUE;
	double[][] lambdaVectors;

	/**
	 * Stores the velocity of each particle
	 */
	private double[][] velocity;
	/**
	 * Stores the personal best solutions found so far for each particle
	 */
//	private Solution pbest_[];// _[];

	/**
	 * Stores the none dominated leaders
	 */
	private CrowdingArchive archive;

	int H_;

	Solution[] indArray_;

	// select the aggregation function to be used
	String functionType_;

	// store the number of the particles' evaluations
	int iteration;
	Operator cloneoperator;
	Operator mutationOperator;
	Operator crossoverOperator;

	int maxIterations;

	private Distance distance_;

	public AgmopsowithR2oldversion(Problem problem) {
		super(problem);
		this.problem = problem;
	} // MOPSOD

	public SolutionSet execute() throws JMException, ClassNotFoundException {

		// to make the algo faster use archiveSize param instead of 100000, this
		// is used here to retrieve as much as possible non-dominated solutions


		functionType_ = "_NBI";
		maxIterations = ((Integer) this.getInputParameter("maxIterations"))
				.intValue();
		populationSize = ((Integer) this.getInputParameter("swarmSize"))
				.intValue();

		archive = new CrowdingArchive(populationSize, problem.getNumberOfObjectives());
		int clonesize = (int) populationSize / 5;

		SolutionSet clonepopulation = new SolutionSet(clonesize);
		int evelations = 0;
		int max_evelations = populationSize * maxIterations;

		iteration = 0;

		population = new SolutionSet(populationSize);
		cpopulation = new SolutionSet(populationSize);
		temppopulation = new SolutionSet(populationSize * 2);
		indArray_ = new Solution[problem.getNumberOfObjectives()];

		// dominance_ = new DominanceComparator();
		distance_ = new Distance();

		cloneoperator = operators_.get("clone");
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");

		H_ = 13; // 23 for 300 and 33 for 595 to be used with 3 objective
		// problems

		T_ = (int) populationSize / 5;
		neighborhood_ = new int[populationSize][T_];
		velocity = new double[this.populationSize][problem
				.getNumberOfVariables()];

		idealPoint = new double[problem.getNumberOfObjectives()];
		nadirPoint = new double[problem.getNumberOfObjectives()];
		lambdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];

		leader_ind = new SolutionSet(populationSize);
		lambdaVectors = new createWeight(problem, populationSize, lambdaVectors).initUniformWeightnorm();
		lambdaVectors0 = deepCopy.deepCopysDouble2d(lambdaVectors);
		initNeighborhood();

		// initialize population
		population = initPopulation();
		population.printObjectivesToFile("AgMOPSO11_variable" + problem.getNumberOfObjectives());
		population.printVariablesToFile("AgMOPSO11_variable_30" + problem.getNumberOfObjectives());
		// initialize the Ideal Point
		initIdealPoint(population);

		//this.createResultsFolders();

		// initialize velocity
		this.initVelocity();
		double[] R2incdicator;
		int[] index;
		index = new int[200];
		SolutionSet temp;
		temp = new SolutionSet(populationSize);
		R2incdicator = new double[temppopulation.size() + archive.size()];
		// STEP 2. Update
		while (evelations < max_evelations) {

			//1.CLONE POPULATION
			distance_.crowdingDistanceAssignment(archive, problem.getNumberOfObjectives());
			archive.sort(new jmetal.util.comparators.CrowdingComparator());
			//get the clone population from the first front
			clonepopulation.clear();
			for (int k = 0; k < archive.size() && k < clonesize; k++) {
				clonepopulation.add(archive.get(k));
			} // for
			cpopulation = (SolutionSet) cloneoperator.execute(clonepopulation);

			temppopulation.clear();
			for (int i = 0; i < cpopulation.size(); i++) {
				Solution[] particle2 = new Solution[2];
				int ran;
				particle2[0] = cpopulation.get(i);
				ran = PseudoRandom.randInt(0, cpopulation.size() - 1);
				particle2[1] = cpopulation.get(ran);

				Solution[] offSpring = (Solution[]) crossoverOperator.execute(particle2);
				mutationOperator.execute(offSpring[0]);
				problem.evaluate(offSpring[0]);
				if (problem.getNumberOfConstraints() != 0)
					problem.evaluateConstraints(offSpring[0]);
				updateReference(offSpring[0]);
				//offSpring[0].setsearch_type(1);
				temppopulation.add(offSpring[0]);
				evelations++;
			}
//            for(int i=0;i<temppopulation.size();i++){
//                archive.add(temppopulation.get(i));
//            }
			for (int i = 0; i < archive.size(); i++) {
				temppopulation.add(archive.get(i));
			}
			R2incdicator = R2__(temppopulation);
			index = sortR2(R2incdicator);
			//get the clone population from the first front

			temp.clear();
			for (int i = 0; i < index.length; i++) {
				if ((R2incdicator[index[i]] != 0) && (i < this.populationSize)) {
					temp.add(temppopulation.get(index[i]));
				} else
					break;
			}
			archive.clear();
			for (int i = 0; i < temp.size(); i++)
				archive.add(temp.get(i));
			iteration++;
			weightVectorAdaption();


			//PSO
			find_leader();
			double speed[][] = this.computeSpeed();
			this.evaluatePopulation(speed);

			temppopulation.clear();
			for (int i = 0; i < archive.size(); i++) {
				temppopulation.add(archive.get(i));
			}
			for (int i = 0; i < population.size(); i++) {
				temppopulation.add(population.get(i));
			}

			R2incdicator = R2__(temppopulation);
			index = sortR2(R2incdicator);
			//get the clone population from the first front

			temp.clear();
			for (int i = 0; i < index.length; i++) {
				if ((R2incdicator[index[i]] != 0) && (i < this.populationSize)) {
					temp.add(temppopulation.get(index[i]));
				} else
					break;
			}
			archive.clear();
			for (int i = 0; i < temp.size(); i++)
				archive.add(temp.get(i));

			evelations += populationSize;
			iteration++;
			weightVectorAdaption();
		}
		return archive;
	}

	/**
	 *
	 */
	public int[] sortR2(double[] R2) {
		int[] index;
		index = new int[R2.length];
		for (int i = 0; i < R2.length; i++) {
			index[i] = i;
		}
		for (int i = 0; i < R2.length; i++) {
			for (int j = 0; j < R2.length - 1; j++) {
				if (R2[index[j]] < R2[index[j + 1]]) {
					int temp = index[j];
					index[j] = index[j + 1];
					index[j + 1] = temp;
				}
			}
		}
		return index;
	}
	public double[] R2__(SolutionSet archive) {
		double[][] TCH;
		TCH = new double[lambdaVectors.length][archive.size()];
		double[] R2indicator;
		R2indicator = new double[archive.size()];
		for (int k = 0; k < lambdaVectors.length; k++) {
			for (int j = 0; j < archive.size(); j++) {
				TCH[k][j] = asf(lambdaVectors[k], archive.get(j));
			}
		}
		int[] index;
		index = new int[lambdaVectors.length];
		double min_;
		for (int i = 0; i < lambdaVectors.length; i++) {
			min_ = 1E+10;
			for (int j = 0; j < archive.size(); j++) {
				if (TCH[i][j] < min_) {
					min_ = TCH[i][j];
					index[i] = j;
				}
			}
			if (min_ == 0) {
				R2indicator[index[i]] = min_ + 0.000001;
			} else {
				R2indicator[index[i]] = min_;
			}
		}

		return R2indicator;
	}

	public double asf(double[] lambdaV, Solution p) {
		double[] temp;
		temp = new double[p.getNumberOfObjectives()];
		for (int i = 0; i < p.getNumberOfObjectives(); i++) {
			temp[i] = Math.abs(p.getObjective(i) - this.idealPoint[i]) / lambdaV[i];
		}
		return Arrays.stream(temp).max().getAsDouble();
	}

	public void find_leader() {
		int best_ind;
		double minFit, fitnesse;
		for (int i = 0; i < this.populationSize; i++) {
			best_ind = -1;
			minFit = Double.MAX_VALUE;
			for (int j = 0; j < archive.size(); j++) {
				fitnesse = this.fitnessFunction(archive.get(j), this.lambdaVectors[i]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					best_ind = j;
				}
			}
			if (fitnessFunction(leader_ind.get(i), lambdaVectors[i]) > minFit) {
				leader_ind.replace(i, new Solution(archive.get(best_ind)));
			}
		}
	}

	public void evaluatePopulation(double[][] speed) throws JMException {

		SolutionSet pop = this.population;

		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle = ;
			for (int var = 0; var < pop.get(n).getDecisionVariables().length; var++) {
				pop.get(n).getDecisionVariables()[var].setValue(pop.get(n).getDecisionVariables()[var].getValue() + speed[n][var]);
				if (pop.get(n).getDecisionVariables()[var].getValue() < problem.getLowerLimit(var)) {
					pop.get(n).getDecisionVariables()[var].setValue(problem.getLowerLimit(var));
					speed[n][var] = 0.0;//speed[n][var] * -1.0;
				}
				if (pop.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
					pop.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
					speed[n][var] = 0.0;//speed[n][var] * -1.0;
				}
			}
		}
		for (int i = 0; i < this.populationSize; i++) {
			Solution particle = pop.get(i);
			// evaluate the new version of the population and update only the
			// particles with better fitness
			problem.evaluate(particle);
			if (problem.getNumberOfConstraints() != 0)
				problem.evaluateConstraints(particle);
			// Update the ideal point
			updateReference(particle);
			// Update of solutions
			updateProblem(particle, i, speed[i]);
			//this.leaders_.add(particle);
//            this.archive.add(particle);
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
				if (j < T_) {
					neighborhood_[k][j] = index.get(j);
				}
			}
		}
	}

	public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
		// list : the set of the indexes of selected mating parents
		// cid  : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
			} else {
				p = PseudoRandom.randInt(0, populationSize - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			//if (flag) list.push_back(p);
			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection

	public boolean updateProblem(Solution indiv, int id, double speed[]) {

		population.replace(id, new Solution(indiv)); // change position
		this.velocity[id] = speed; // update speed

		//
		return true;

	} // updateProblem
	// /////////////////////////////////////////////////////

	private double[][] computeSpeed() throws JMException {
		double W;
		double[][] speed = new double[this.populationSize][problem.getNumberOfVariables()];

		int l2;

		Variable[] pbest, lbest, gbest;
		for (int n = 0; n < this.population.size(); n++) {
			Variable[] particle = population.get(n).getDecisionVariables();
			double f;

			pbest = leader_ind.get(n).getDecisionVariables();

			l2 = PseudoRandom.randInt(0, this.archive.size() - 1); // select random  leader
			gbest = archive.get(l2).getDecisionVariables();

			Vector<Integer> p = new Vector<Integer>();
			matingSelection(p, n, 1, 1); //Select a neighborhood sub-problem of n
			lbest = leader_ind.get(p.get(0)).getDecisionVariables();
			f = 0.5;//diversity(leader_ind.get(n),lambdaVectors[n])/max_d;
			double c = PseudoRandom.randDouble();//diversity(leader_ind.get(n),lambdaVectors[n])/max_d;
			W = PseudoRandom.randDouble(0.1, 0.5);
			for (int var = 0; var < particle.length; var++) {
				speed[n][var] = (W * velocity[n][var])
						+ c * (pbest[var].getValue() - particle[var].getValue())
						+ f * (lbest[var].getValue() - gbest[var].getValue());
			}// end for
		}
		return speed;
	}

	// ///////////////////////////////////////////////////////////////////////
	//��ʼ�������ٶ�
	private void initVelocity() {
		for (int i = 0; i < this.populationSize; i++)
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				velocity[i][j] = 0.0;
			}
	}

	// ////////////////////////////////////////////////////////////////////////
	//��ʼ����Ⱥ��α������
	public SolutionSet initPopulation() throws JMException,
			ClassNotFoundException {
		SolutionSet pop = new SolutionSet(this.populationSize);
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem);
			problem.evaluate(newSolution);
			if (this.problem.getNumberOfConstraints() != 0)
				problem.evaluateConstraints(newSolution);
			// evaluations++;
			pop.add(newSolution);
			leader_ind.add(newSolution);
			archive.add(newSolution);
		}
		return pop;
	} // initPopulation
	// ///////////////////////////////////////////////////////////////////////////

	// ******************************************************************
	//��ʼ���ο���
	void initIdealPoint(SolutionSet pop) throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = 1.0e+30;
			indArray_[i] = new Solution(problem);
			problem.evaluate(indArray_[i]);
			// evaluations++;
		} // for

		for (int i = 0; i < populationSize; i++)
			updateReference(pop.get(i));//����Ⱥ�нϺõĸ���������ǰ��������

	} // initIdealPoint

	// *******************************************************************
	//���²ο���
	double fitnessFunction(Solution indiv, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_NBI")) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

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
		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation

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
	void updateReference(Solution individual) {
		for (int n = 0; n < problem.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < idealPoint[n]) {
				idealPoint[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference
} // MOPSOD