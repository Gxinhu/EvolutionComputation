package jmetal.metaheuristics.agmopso;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.createWeight;
import jmetal.util.savesomething.savetofile;

import java.util.Vector;

public class AgMOPSO extends Algorithm {

	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private QualityIndicator indicator;
	int run;
	int t;
	int[][] neighborhood;
	public String curDir = System.getProperty("user.dir");
	private double[][] realtimeIGD;
	private double[][] realtimeSpeard;
	private double[][] realtimeGD;
	private int populationSize;
	/**
	 * Stores the population
	 */
	private SolutionSet population, temppopulation, cpopulation, leader_ind;
	/**
	 * Z vector (ideal point)
	 */
	double[] idealPoint;
	/**
	 * Lambda vectors
	 */
	private double max_d = Double.MIN_VALUE;
	double[][] lamdaVectors;

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

	int h;

	private Solution[] indarray;

	// select the aggregation function to be used
	private String functiontype;

	// store the number of the particles' evaluations
	int iteration;
	Operator cloneoperator;
	Operator mutationOperator;
	Operator crossoverOperator;
	int maxIterations;
	private Distance distance;
	private int runtimes;
	private boolean save;

	public AgMOPSO(Problem problem, QualityIndicator indicator, int i, boolean save) {
		super(problem);
		this.problem = problem;
		this.indicator = indicator;
		this.runtimes = i;
		this.save = save;
	} // MOPSOD

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		// to make the algo faster use archiveSize param instead of 100000, this
		// is used here to retrieve as much as possible non-dominated solutions

		functiontype = "_NBI";
		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = ((Integer) this.getInputParameter("swarmSize"))
				.intValue();
		realtimeIGD = new double[maxIterations / 10 + 1][2];
		realtimeSpeard = new double[maxIterations / 10 + 1][2];
		realtimeGD = new double[maxIterations / 10 + 1][2];
		archive = new CrowdingArchive(populationSize, problem.getNumberOfObjectives());
		int clonesize = populationSize / 5;

		SolutionSet clonepopulation = new SolutionSet(clonesize);
		int evelations = 0;
		int max_evelations = populationSize * maxIterations;

		iteration = 0;

		population = new SolutionSet(populationSize);
		cpopulation = new SolutionSet(populationSize);
		temppopulation = new SolutionSet(populationSize);
		indarray = new Solution[problem.getNumberOfObjectives()];

		// dominance_ = new DominanceComparator();
		distance = new Distance();

		cloneoperator = operators_.get("clone");
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");

		t = populationSize / 5;
		neighborhood = new int[populationSize][t];
		velocity = new double[this.populationSize][problem
				.getNumberOfVariables()];

		idealPoint = new double[problem.getNumberOfObjectives()];

		lamdaVectors = new double[populationSize][problem
				.getNumberOfObjectives()];

		leader_ind = new SolutionSet(populationSize);
		lamdaVectors = new createWeight(problem, populationSize, lamdaVectors).initUniformWeightWs();
		initNeighborhood();

		// initialize population
		population = initPopulation();
		// initialize the Ideal Point
		initIdealPoint(population);
		calulateindicator();
		// initialize velocity
		this.initVelocity();
		// STEP 2. Update
		while (iteration <= maxIterations) {
			problem_.dynamicChange(iteration);
			for (int i = 0; i < population.size(); i++) {
				double temp = population.get(i).getObjective(1);
				problem.evaluate(population.get(i));
				double temp1 = population.get(i).getObjective(1);
				if (temp != temp1) {
					System.out.println(iteration);
					break;
				}
			}
			//1.CLONE POPULATION
			distance.crowdingDistanceAssignment(archive, problem.getNumberOfObjectives());
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
				if (problem.getNumberOfConstraints() != 0) {
					problem.evaluateConstraints(offSpring[0]);
				}
				updateReference(offSpring[0]);
				temppopulation.add(offSpring[0]);
				evelations++;
			}
			for (int i = 0; i < temppopulation.size(); i++) {
				archive.add(temppopulation.get(i));
			}
			iteration++;
			problem_.dynamicChange(iteration);
			if (iteration >= maxIterations) {
				break;
			}
			find_leader();

			double[][] speed = this.computeSpeed();
			this.evaluatePopulation(speed);
			evelations += populationSize;
			iteration++;
			if (iteration % 10 == 0) {
				calulateindicator();
			}

		}

//		写入csv文件,随着迭代次数的指标变化
		if (save) {
			savetofile savetofile = new savetofile(problem, "out/IGD/" + problem.getName(), runtimes, realtimeIGD);
			savetofile.save();
			savetofile = new savetofile(problem, "out/GD/" + problem.getName(), runtimes, realtimeGD);
			savetofile.save();
			savetofile = new savetofile(problem, "out/Spread/" + problem.getName(), runtimes, realtimeSpeard);
			savetofile.save();
		}
		return archive;
	}

	/**
	 * caculate the indicator with time
	 */
	private void calulateindicator() {
		if (this.save) {
			realtimeSpeard[iteration / 10][0] = iteration;
			realtimeSpeard[iteration / 10][1] = indicator.getGeneralizedSpread(archive);
			realtimeIGD[iteration / 10][0] = iteration;
			realtimeIGD[iteration / 10][1] = indicator.getCEC_IGD(archive);
			realtimeGD[iteration / 10][0] = iteration;
			realtimeGD[iteration / 10][1] = indicator.getGD(archive);
		}
	}

	/**
	 *
	 */

	public void find_leader() {
		int bestInd;
		double minFit, fitnesse;
		for (int i = 0; i < this.populationSize; i++) {
			bestInd = -1;
			minFit = Double.MAX_VALUE;
			for (int j = 0; j < archive.size(); j++) {
				fitnesse = this.fitnessFunction(archive.get(j), this.lamdaVectors[i]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					bestInd = j;
				}
			}
			if (fitnessFunction(leader_ind.get(i), lamdaVectors[i]) > minFit) {
				leader_ind.replace(i, new Solution(archive.get(bestInd)));
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
					speed[n][var] = 0.0;
					//speed[n][var] * -1.0;
				}
				if (pop.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
					pop.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
					speed[n][var] = 0.0;
					//speed[n][var] * -1.0;
				}
			}
		}
		for (int i = 0; i < this.populationSize; i++) {
			Solution particle = pop.get(i);
			// evaluate the new version of the population and update only the
			// particles with better fitness
			problem.evaluate(particle);
			if (problem.getNumberOfConstraints() != 0) {
				problem.evaluateConstraints(particle);
			}
			// Update the ideal point
			updateReference(particle);
			// Update of solutions
			updateProblem(particle, i, speed[i]);
			//this.leaders_.add(particle);
			this.archive.add(particle);
		}

	}

	public void initNeighborhood() {
		double[] x = new double[populationSize];
		int[] idx = new int[populationSize];

		for (int i = 0; i < populationSize; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize; j++) {
				x[j] = Utils.distVector(lamdaVectors[i], lamdaVectors[j]);
				//x[j] = dist_vector(population[i].namda,population[j].namda);
				idx[j] = j;
				//System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]: "+idx[j]) ;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize, t);
			//minfastsort(x,idx,population.size(),niche);

			System.arraycopy(idx, 0, neighborhood[i], 0, t);
		} // for
	} // initNeighborhood

	public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
		// list : the set of the indexes of selected mating parents
		// cid  : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood[cid][r];
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

	public boolean updateProblem(Solution indiv, int id, double[] speed) {

		population.replace(id, new Solution(indiv)); // change position
		this.velocity[id] = speed; // update speed

		//
		return true;

	} // updateProblem
	// /////////////////////////////////////////////////////

	private double[][] computeSpeed() throws JMException {
		double w;
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
			w = PseudoRandom.randDouble(0.1, 0.5);
			for (int var = 0; var < particle.length; var++) {
				speed[n][var] = (w * velocity[n][var])
						+ c * (pbest[var].getValue() - particle[var].getValue())
						+ f * (lbest[var].getValue() - gbest[var].getValue());
//						+ f*( gbest[var].getValue()-particle[var].getValue());			// end for
			}
		}

		return speed;
	}

	// ///////////////////////////////////////////////////////////////////////

	private void initVelocity() {
		for (int i = 0; i < this.populationSize; i++) {
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				velocity[i][j] = 0.0;
			}
		}
	}

	// ////////////////////////////////////////////////////////////////////////

	public SolutionSet initPopulation() throws JMException,
			ClassNotFoundException {
		SolutionSet pop = new SolutionSet(this.populationSize);
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem);
			problem.evaluate(newSolution);
			if (this.problem.getNumberOfConstraints() != 0) {
				problem.evaluateConstraints(newSolution);
			}
			// evaluations++;
			pop.add(newSolution);
			leader_ind.add(newSolution);
			archive.add(newSolution);
		}
		return pop;
	} // initPopulation
	// ///////////////////////////////////////////////////////////////////////////

	// ******************************************************************

	private void initIdealPoint(SolutionSet pop) throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = 1.0e+30;
			indarray[i] = new Solution(problem);
			problem.evaluate(indarray[i]);
			// evaluations++;
		} // for

		for (int i = 0; i < populationSize; i++) {
			updateReference(pop.get(i));//����Ⱥ�нϺõĸ���������ǰ��������
		}

	} // initIdealPoint


	private double fitnessFunction(Solution indiv, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functiontype.equals("_NBI")) {
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
					+ functiontype);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation

	void updateReference(Solution individual) {
		for (int n = 0; n < problem.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < idealPoint[n]) {
				idealPoint[n] = individual.getObjective(n);

				indarray[n] = individual;
			}
		}
	} // updateReference
} // MOPSOD