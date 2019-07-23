package jmetal.metaheuristics.mmopso;

import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.CrowdingComparator;

import java.io.File;

public class MMOPSO extends Algorithm {

	public MMOPSO(Problem problem) {
		//this.run = r;
		//this.problem = problem;
		super(problem);
		// TODO Auto-generated constructor stub
	}

	private static final long serialVersionUID = 2107684627645440737L;
	//private Problem problem;

	int run;
	public String curDir = System.getProperty("user.dir");

	private int populationSize;
	/**
	 * Stores the population
	 */
	private SolutionSet population, temppopulation;
	/**
	 * Z vector (ideal point)
	 */
	double[] idealPoint;
	/**
	 * Lambda vectors
	 */

	double[][] lamdaVectors;

	int[] leader_ind;

	/**
	 * Stores the maximum size for the archive
	 */
	private int selectionArchiveSize;
	/**
	 * Stores the velocity of each particle
	 */
	private double[][] velocity;
	/**
	 * Stores the personal best solutions found so far for each particle
	 */
	private Solution[] pbest_;// _[];

	/**
	 * Stores the none dominated leaders
	 */
	private CrowdingArchive archive;
	//private CrowdingArchiveBasedonRanking leaders_;
	int H_;
	int neighbourhoodSize = 20;
	Solution[] indArray_;

	// select the aggregation function to be used
	String functionType_;

	// store the number of the particles' evaluations
	int iteration;
	Operator mutationOperator;
	Operator crossoverOperator;

	int maxIterations;

	private Distance distance_;
	int nr_;

	/*public MMOPSO(Problem problem, int r) {
		this.run = r;
		this.problem = problem;

	} // MOPSOD*/

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		// to make the algo faster use archiveSize param instead of 100000, this
		// is used here to retrieve as much as possible non-dominated solutions


		functionType_ = ((String) getInputParameter("functionType"));
		selectionArchiveSize = ((Integer) getInputParameter("selectionArchiveSize"))
				.intValue();
		maxIterations = ((Integer) this.getInputParameter("maxIterations"))
				.intValue();
		populationSize = ((Integer) this.getInputParameter("swarmSize"))
				.intValue();
		archive = new CrowdingArchive(populationSize, problem_.getNumberOfObjectives());
		int evelations = 0;
		int max_evelations = populationSize * maxIterations;
		double probability = 0;
		int mutated_num;
		int good_mutated2;

		iteration = 0;
		pbest_ = new Solution[this.populationSize];
		population = new SolutionSet(populationSize);
		temppopulation = new SolutionSet(populationSize);
		indArray_ = new Solution[problem_.getNumberOfObjectives()];
		//leaders_ = new CrowdingArchiveBasedonRanking(selectionArchiveSize,
		//		problem.getNumberOfObjectives(), problem.getNumberOfVariables());

		// dominance_ = new DominanceComparator();
		distance_ = new Distance();

		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");

		H_ = 33; // 23 for 300 and 33 for 595 to be used with 3 objective
		// problems

		velocity = new double[this.populationSize][problem_
				.getNumberOfVariables()];

		idealPoint = new double[problem_.getNumberOfObjectives()];

		lamdaVectors = new double[populationSize][problem_
				.getNumberOfObjectives()];

		leader_ind = new int[populationSize];
		initUniformWeight();

		// initialize population
		population = initPopulation();

		// initialize the Ideal Point
		initIdealPoint(population);

		//this.createResultsFolders();

		// initialize velocity
		this.initVelocity();

		// assign lambda vector to each particle
		this.orderPopulation(population);

		// STEP 2. Update
		probability = 0.9;
		good_mutated2 = 0;
		mutated_num = 0;
		while (evelations < max_evelations) {
			find_leader();
			double[][] speed = this.computeSpeed(probability);
			this.evaluatePopulation(speed);

			distance_.crowdingDistanceAssignment(archive,
					problem_.getNumberOfObjectives());
			archive.sort(new CrowdingComparator());
			iteration++;
			//System.out.println("ev"+evelations+"archive size"+archive.size());
			mutated_num = archive.size();
			good_mutated2 = archive.size() / 4;
			temppopulation.clear();
			//mutated_num=0;
			for (int i = 0; i < mutated_num; i++) {
				Solution[] particle2 = new Solution[2];
				int ran;
				particle2[0] = archive.get(i);
				ran = PseudoRandom.randInt(0, good_mutated2 - 1);
				particle2[1] = archive.get(ran);
				Solution[] offSpring = (Solution[]) crossoverOperator.execute(particle2);
				mutationOperator.execute(offSpring[0]);
				problem_.evaluate(offSpring[0]);
				if (problem_.getNumberOfConstraints() != 0) {
					problem_.evaluateConstraints(offSpring[0]);
				}
				updateReference(offSpring[0]);
				//offSpring[0].setsearch_type(1);
				temppopulation.add(offSpring[0]);
				evelations++;
			}
			for (int i = 0; i < temppopulation.size(); i++) {
				archive.add(temppopulation.get(i));
			}
			//System.out.println("ev:"+evelations+" mutated num:"+mutated_num+" good_mutated:"+good_mutated);
			//evelations+=archive.size();
			//}
			evelations += populationSize;

			// leaders_.printObjectivesToFile(curDir + "/RESULTS-DDMOPSO-" +
			// problem.getName() + "/RUN-" + run + "/FUN/Leaders.GEN-" +
			// iteration); leaders_.printVariablesToFile(curDir +
			// "/RESULTS-DDMOPSO-" + problem.getName() + "/RUN-" + run +
			// "/VAR/Leaders.GEN-" + iteration);
			//
			// this.population.printObjectivesToFile(curDir +
			// "/RESULTS-DDMOPSO-" + problem.getName() + "/RUN-" + run +
			// "/FUN/Swarm.GEN-" + iteration);
			// this.population.printVariablesToFile(curDir + "/RESULTS-DDMOPSO-"
			// + problem.getName() + "/RUN-" + run + "/" + "/VAR/Swarm.GEN-" +
			// iteration);

			//  this.archive.printObjectivesToFile(curDir + "/RESULTS-DDMOPSO-" +
			//  problem.getName() + "/RUN-" + run + "/FUN/eArchive.GEN-" +
			//  iteration); 
			// this.archive.printVariablesToFile(curDir +
			// "/RESULTS-DDMOPSO-" + problem.getName() + "/RUN-" + run + "/" +
			// "/VAR/eArchive.GEN-" + iteration);
			// System.out.println(iteration);
		}
		//archive.printObjectivesToFile(curDir + "/RESULTS-DDMOPSO-"
		//		+ problem.getName() + "/RUN-" + run + "-FUN");
		//archive.printVariablesToFile(curDir + "/RESULTS-DDMOPSO-"
		//		+ problem.getName() + "/RUN-" + run + "-VAR");
		//System.out.println("ev:"+evelations+" total Search 1:"+good_sum2[0]+" total Search 2:"+good_sum2[1]);
		//System.out.println("ev:"+evelations+" mutated num:"+mutated_num2+" good_mutated:"+good_mutated2);

		return archive;
	}

	/**
	 *
	 */

	public void find_leader() {
		int best_ind;
		double minFit, fitnesse;
		for (int i = 0; i < this.populationSize; i++) {
			best_ind = -1;
			minFit = Double.MAX_VALUE;
			for (int j = 0; j < archive.size(); j++) {
				fitnesse = this.fitnessFunction(archive.get(j), this.lamdaVectors[i]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					best_ind = j;
				}
			}
			leader_ind[i] = best_ind;
		}
	}

	public void orderPopulation(SolutionSet pop) {
		population = new SolutionSet(populationSize);

		double[][] fitnesses = new double[this.populationSize][this.populationSize];
		for (int i = 0; i < this.populationSize; i++) {
			for (int j = 0; j < this.populationSize; j++) {
				fitnesses[i][j] = this.fitnessFunction(pop.get(i),
						this.lamdaVectors[j]);
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
			this.pbest_[i] = new Solution(pop.get(particleIndex));
			//this.leaders_.add(pop.get(particleIndex));
			this.archive.add(pop.get(particleIndex));
		}

	}

	/**
	 * Update the position of each particle
	 *
	 * @throws JMException
	 */
	private SolutionSet computeNewPositions(double[][] speed)
			throws JMException {

		SolutionSet pop = this.population;

		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle = ;
			for (int var = 0; var < pop.get(n).getDecisionVariables().length; var++) {
				pop.get(n).getDecisionVariables()[var]
						.setValue(pop.get(n).getDecisionVariables()[var]
								.getValue() + speed[n][var]);
				if (pop.get(n).getDecisionVariables()[var].getValue() < problem_
						.getLowerLimit(var)) {
					pop.get(n).getDecisionVariables()[var].setValue(problem_
							.getLowerLimit(var));
					speed[n][var] = speed[n][var] * -1.0;
				}
				if (pop.get(n).getDecisionVariables()[var].getValue() > problem_
						.getUpperLimit(var)) {
					pop.get(n).getDecisionVariables()[var].setValue(problem_
							.getUpperLimit(var));
					speed[n][var] = speed[n][var] * -1.0;
				}
			}
		}
		return pop;
	} // computeNewPositions
	// //////////////////////////////////////////////////////////////////////////////////////

	public void evaluatePopulation(double[][] speed) throws JMException {

		SolutionSet pop = this.computeNewPositions(speed);
		for (int i = 0; i < this.populationSize; i++) {
			Solution particle = pop.get(i);
			// evaluate the new version of the population and update only the
			// particles with better fitness
			problem_.evaluate(particle);
			if (problem_.getNumberOfConstraints() != 0) {
				problem_.evaluateConstraints(particle);
			}
			// Update the ideal point
			updateReference(particle);
			// Update of solutions
			updateProblem(particle, i, speed[i]);
			//this.leaders_.add(particle);
			this.archive.add(particle);
		}

	}

	// ///////////////////////////////////////////////////////////////////////////////////////
	public void initUniformWeight() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			for (int n = 0; n < populationSize; n++) {
				double a = 1.0 * n / (populationSize - 1);
				lamdaVectors[n][0] = a;
				lamdaVectors[n][1] = 1 - a;
				nw++;
			} // for
		} // if
		else {
			int i, j;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						lamdaVectors[nw][0] = (1.0 * i) / H_;
						lamdaVectors[nw][1] = (1.0 * j) / H_;
						lamdaVectors[nw][2] = 1.0 * (H_ - i - j) / H_;
						nw++;
					} // if
				} // for
			} // for
		} // else

		if (nw != populationSize) {
			System.out.println(nw + "---" + (populationSize));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
	} // initUniformWeight
	// ////////////////////////////////////////////////////////////////////////////////////////

	public boolean updateProblem(Solution indiv, int id, double[] speed) {

		population.replace(id, new Solution(indiv)); // change position
		this.velocity[id] = speed; // update speed

		//
		return true;

	} // updateProblem
	// /////////////////////////////////////////////////////

	private double[][] computeSpeed(double pro) throws JMException {
		double r2, W, C2;
		double[][] speed = new double[this.populationSize][problem_
				.getNumberOfVariables()];
		int l2;

		Variable[] gbest;
		for (int n = 0; n < this.population.size(); n++) {

			// nbest = population.get(l1).getDecisionVariables();//
			// pbest_[p].getDecisionVariables();
			Variable[] particle = population.get(n).getDecisionVariables();

			double rand = PseudoRandom.randDouble();
			if (rand < pro) {
				//l2 = this.leaderSelection2(n); // select leader based on
				// decomposition
				//l2 = PseudoRandom.randInt(0, neighbourhoodSize-1);
				//l2 = leader_ind[neighborhood_[n][l2]]; // select leader based on
				// decomposition
				l2 = leader_ind[n];
			} else {
				l2 = PseudoRandom.randInt(0, this.archive.size() - 1); // select
				// random
				// leader
			}

			gbest = archive.get(l2).getDecisionVariables();

			for (int var = 0; var < particle.length; var++) {

				//r1 = PseudoRandom.randDouble();
				r2 = PseudoRandom.randDouble();
				//C1 = PseudoRandom.randDouble(1.5, 2.0);
				C2 = PseudoRandom.randDouble(1.5, 2.0);
				W = PseudoRandom.randDouble(0.1, 0.5);

				//Variable[] pbest = pbest_[n].getDecisionVariables();

				//rand = PseudoRandom.randDouble();
				speed[n][var] = (W * velocity[n][var]) +
						+C2 * r2 * (gbest[var].getValue() - particle[var].getValue());
						/*speed[n][var] = (W * velocity[n][var]) +
						+ C2 * r2
						* (gbest[var].getValue() - particle[var].getValue());*/

			}// end for

		}

		return speed;
	}

	// ///////////////////////////////////////////////////////////////////////
	private void initVelocity() {
		for (int i = 0; i < this.populationSize; i++) {
			for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
				velocity[i][j] = 0.0;
			}
		}
	}

	// ////////////////////////////////////////////////////////////////////////
	public SolutionSet initPopulation() throws JMException,
			ClassNotFoundException {
		SolutionSet pop = new SolutionSet(this.populationSize);
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			if (this.problem_.getNumberOfConstraints() != 0) {
				problem_.evaluateConstraints(newSolution);
			}
			// evaluations++;
			pop.add(newSolution);
		}
		return pop;
	} // initPopulation
	// ///////////////////////////////////////////////////////////////////////////


	public int leaderSelection2(int cid) {

		int p = 0;

		p = PseudoRandom.randInt(0, archive.size() - 1);

		int leaderIndex = 0;
		double fit = Double.MAX_VALUE;
		for (int i = 0; i < archive.size(); i++) {
			Solution individua = this.archive.get(i);
			double f = this.fitnessFunction(individua, lamdaVectors[cid]);
			if (f < fit) {
				fit = f;
				leaderIndex = i;
			}
		}
		p = leaderIndex;

		return p;

	} // leaders selection

	// ******************************************************************
	void initIdealPoint(SolutionSet pop) throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			idealPoint[i] = 1.0e+30;
			indArray_[i] = new Solution(problem_);
			problem_.evaluate(indArray_[i]);
			// evaluations++;
		} // for

		for (int i = 0; i < populationSize; i++) {
			updateReference(pop.get(i));
		}

	} // initIdealPoint

	// ***************************************************************

	double fitnessFunction(Solution individual, double[] lamda) {

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				// double diff = Math.abs(individual.getObjective(n)
				// - this.idealPoint[n]);

				double diff = Math.abs(individual.getObjective(n)
						- idealPoint[n]);

				double feval;
				if (lamda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lamda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			return maxFun;
		} // if

		else if (functionType_.equals("_WSUM")) {

			double sum = 0;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				sum += (lamda[n]) * individual.getObjective(n);
			}
			return sum;

		} // if
		else if (functionType_.equals("_NBI")) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

			d1 = d2 = nl = 0.0;
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - idealPoint[i]) * lamda[i];
				nl += (lamda[i] * lamda[i]);
			}
			if (d1 > 0) {
				d1 = d1 / Math.sqrt(nl);
			} else {
				d1 = -d1 / Math.sqrt(nl);
			}
			//d1 = Math.abs(d1) / Math.sqrt(nl);
			/*if (nl == 0.0) {
				System.out
						.println("ERROR: dived by zero(bad weihgted vector)\n");
				System.exit(0);
			}*/
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				//d2 += Math.pow((individual.getObjective(i) - idealPoint[i])
				//		- (d1 * lamda[i]), 2.0);
				fin = (individual.getObjective(i) - idealPoint[i]) - (d1 * lamda[i]);
				//d2 += (individual.getObjective(i) - idealPoint[i]) - (d1 * lamda[i]), 2.0);
				d2 += (fin * fin);
			}
			d2 = Math.sqrt(d2);
			fin = (d1 + theta * d2);
			return fin;
		} else {
			System.out.println("SDMOPSO.fitnessFunction: unknown type "
					+ functionType_);
			return 0;
		}

	} // fitnessEvaluation
	/*double fitnessFunction(Solution individual, double[] lamda) {

		int i;
		double d1, d2, nl;
		double theta = 5.0;
		double fin;

		d1 = d2 = nl = 0.0;
		for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
			d1 += (individual.getObjective(i) - idealPoint[i]) * lamda[i];
			nl += (lamda[i]*lamda[i]);
		}
		nl = Math.sqrt(nl);
		if(d1>0){
			d1 = d1 / nl;
		}
		else
			d1 = -d1 / nl;
		//d1 = Math.abs(d1) / Math.sqrt(nl);
		if (nl == 0.0) {
			System.out
					.println("ERROR: dived by zero(bad weihgted vector)\n");
			System.exit(0);
		}
		for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
			//d2 += Math.pow((individual.getObjective(i) - idealPoint[i])
			//		- (d1 * lamda[i]), 2.0);
			fin =(individual.getObjective(i) - idealPoint[i]) - (d1 * lamda[i]);
			//d2 += (individual.getObjective(i) - idealPoint[i]) - (d1 * lamda[i]), 2.0);
			d2 += (fin*fin);
		}
		d2 = Math.sqrt(d2);
		fin = (d1 + theta * d2);
		return fin;

	} // fitnessEvaluation*/

	// *******************************************************************
	void updateReference(Solution individual) {
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < idealPoint[n]) {
				idealPoint[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference
	// ///////////////////////////////////////////////////////////////////////

	private void createFolder(String str) {

		boolean success = (new File(str)).mkdir();

		if (success) {
			System.out.println("Directory: " + str + " created");
		} else {

			System.out.println("Directory: " + str + " NOT created");

		}

	}

	// ///////////////////////////////////////////////////////////////////////
	private void createResultsFolders() {

		String str = curDir + "/RESULTS-DDMOPSO-" + problem_.getName();

		str += "/RUN-" + run + "";
		this.createFolder(str);

		String fun = str + "/FUN";
		this.createFolder(fun);

		String var = str + "/VAR";
		this.createFolder(var);

	}

	// ///////////////////////////////////////////////////////////////////////
	private boolean deleteFolder(File dir) {

		if (dir.isDirectory()) {
			String[] children = dir.list();
			for (int i = 0; i < children.length; i++) {
				boolean success = deleteFolder(new File(dir, children[i]));
				if (!success) {
					return false;
				}
			}
		}

		// The directory is now empty so delete it
		return dir.delete();

	}
	// //////////////////////////////////////////////////////////////////

} // MOPSOD