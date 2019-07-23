package jmetal.metaheuristics.smpso;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.wrapper.XReal;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Comparator;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 */
public class CMPSO extends Algorithm {

	/**
	 * Stores the number of particles_ used
	 */
	private int swarmSize_;
	/**
	 * Stores the maximum size for the archive
	 */
	private int all_swarmSize_;
	private int archiveSize_;

	private int function_ev;
	private int max_function_ev;
	private SolutionSet particles_;
	/**
	 * Stores the best_ solutions founds so far for each particles
	 */
	private Solution[] pbest_;

	private Solution[] gbest_;
	/**
	 * Stores the leaders_
	 */

	private SolutionSet archive_;

	/**
	 * Stores the speed_ of each particle
	 */
	private double[][] speed_;
	/**
	 * Stores a comparator for checking dominance
	 */
	private Comparator dominance_;
	/**
	 * Stores a comparator for crowding checking
	 */
	/**
	 * Stores a <code>Distance</code> object
	 */
	private Distance distance_;
	/**
	 * Stores a operator for non uniform mutations
	 */

	QualityIndicator indicators_; // QualityIndicator object

	double ChVel1_;
	double ChVel2_;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public CMPSO(Problem problem) {
		super(problem);

		ChVel1_ = -1;
		ChVel2_ = -1;
	} // Constructor

	public CMPSO(Problem problem, Vector<Double> variables,
	             String trueParetoFront) throws FileNotFoundException {
		super(problem);

		ChVel1_ = -1;
		ChVel2_ = -1;

		hy_ = new Hypervolume();
		jmetal.qualityIndicator.util.MetricsUtil mu = new jmetal.qualityIndicator.util.MetricsUtil();
		trueFront_ = mu.readNonDominatedSolutionSet(trueParetoFront);
		trueHypervolume_ = hy_.hypervolume(
				trueFront_.writeObjectivesToMatrix(),
				trueFront_.writeObjectivesToMatrix(),
				problem_.getNumberOfObjectives());

	} // SMPSO

	private double trueHypervolume_;
	private Hypervolume hy_;
	private SolutionSet trueFront_;
	private double[] deltaMax_;
	private double[] deltaMin_;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public CMPSO(Problem problem, String trueParetoFront)
			throws FileNotFoundException {
		super(problem);
		hy_ = new Hypervolume();
		jmetal.qualityIndicator.util.MetricsUtil mu = new jmetal.qualityIndicator.util.MetricsUtil();
		trueFront_ = mu.readNonDominatedSolutionSet(trueParetoFront);
		trueHypervolume_ = hy_.hypervolume(
				trueFront_.writeObjectivesToMatrix(),
				trueFront_.writeObjectivesToMatrix(),
				problem_.getNumberOfObjectives());

		// Default configuration
		ChVel1_ = -1;
		ChVel2_ = -1;
	} // Constructor

	/**
	 * Initialize all parameter of the algorithm
	 */
	public void initParams() {
		swarmSize_ = ((Integer) getInputParameter("swarmSize")).intValue();
		archiveSize_ = ((Integer) getInputParameter("archiveSize")).intValue();
		max_function_ev = ((Integer) getInputParameter("max_function_ev"))
				.intValue();

		indicators_ = (QualityIndicator) getInputParameter("indicators");
		all_swarmSize_ = swarmSize_ * problem_.getNumberOfObjectives();

		//polynomialMutation_ = operators_.get("mutation");

		function_ev = 0;

		particles_ = new SolutionSet(all_swarmSize_);
		pbest_ = new Solution[all_swarmSize_];
		gbest_ = new Solution[problem_.getNumberOfObjectives()];
		archive_ = new SolutionSet(all_swarmSize_ + 2 * archiveSize_);

		// Create comparators for dominance and crowding distance
		dominance_ = new DominanceComparator();
		distance_ = new Distance();

		// Create the speed_ vector
		speed_ = new double[all_swarmSize_][problem_.getNumberOfVariables()];

		deltaMax_ = new double[problem_.getNumberOfVariables()];
		deltaMin_ = new double[problem_.getNumberOfVariables()];
		for (int i = 0; i < problem_.getNumberOfVariables(); i++) {
			deltaMax_[i] = (problem_.getUpperLimit(i) - problem_
					.getLowerLimit(i)) / 5.0;
			deltaMin_[i] = -deltaMax_[i];
		} // for
	} // initParams

	/**
	 * Update the speed of each particle
	 *
	 * @throws JMException
	 */
	private void computeSpeed(int iter, int miter) throws JMException,
			IOException {
		double r1, r2, r3, W, C1, C2, C3;
		XReal particle, bestParticle, gbestParticle, abestParticle;
		C1 = 4.0 / 3;
		C2 = 4.0 / 3;
		C3 = 4.0 / 3;
		W = 0.9 - 0.5 * function_ev / max_function_ev;
		for (int i = 0; i < particles_.size(); i++) {
			particle = new XReal(particles_.get(i));
			bestParticle = new XReal(pbest_[i]);
			int j = i / swarmSize_;
			gbestParticle = new XReal(gbest_[j]);

			// Select a global best_ from archive for calculating the speed of particle i,
			// bestGlobal
			Solution one;
			int pos1 = PseudoRandom.randInt(0, archive_.size() - 1);
			one = archive_.get(pos1);
			abestParticle = new XReal(one);


			for (int var = 0; var < particle.getNumberOfDecisionVariables(); var++) {
				r1 = PseudoRandom.randDouble(0, 1);
				r2 = PseudoRandom.randDouble(0, 1);
				r3 = PseudoRandom.randDouble(0, 1);
				// Computing the velocity of this particle
				speed_[i][var] = W * speed_[i][var] + C1 * r1 * (bestParticle.getValue(var) - particle.getValue(var)) +
						C2 * r2 * (gbestParticle.getValue(var) - particle.getValue(var)) +
						C3 * r3 * (abestParticle.getValue(var) - particle.getValue(var));
				if (speed_[i][var] > deltaMax_[var]) //constraint the velocity between deltaMax and deltaMin
				{
					speed_[i][var] = deltaMax_[var];
				}
				if (speed_[i][var] < deltaMin_[var]) {
					speed_[i][var] = deltaMin_[var];
				}
			}
		}
	} // computeSpeed

	/**
	 * Update the position of each particle
	 *
	 * @throws JMException
	 */
	private void computeNewPositions() throws JMException {
		for (int i = 0; i < particles_.size(); i++) {
			XReal particle = new XReal(particles_.get(i));
			for (int var = 0; var < particle.getNumberOfDecisionVariables(); var++) {
				particle.setValue(var, particle.getValue(var) + speed_[i][var]);

				if (particle.getValue(var) < problem_.getLowerLimit(var)) {
					particle.setValue(var, problem_.getLowerLimit(var));
					speed_[i][var] = speed_[i][var] * ChVel1_; //
				}
				if (particle.getValue(var) > problem_.getUpperLimit(var)) {
					particle.setValue(var, problem_.getUpperLimit(var));
					speed_[i][var] = speed_[i][var] * ChVel2_; //
				}
			}
		}
	} // computeNewPositions	

	public SolutionSet Getnondominated(SolutionSet solutionSet) {
		SolutionSet ranking_;
		Solution a1, a2;
		int Set_size = solutionSet.size(), num_obj = solutionSet.get(0).getNumberOfObjectives();
		int[] dominateMe = new int[Set_size];

		for (int i = 0; i < Set_size; i++) {
			dominateMe[i] = 0;
		}
		ranking_ = new SolutionSet(solutionSet.getMaxSize());

		// iDominate[k] contains the list of solutions dominated by k
		for (int i = 0; i < Set_size; i++) {
			a1 = solutionSet.get(i);
			if (dominateMe[i] == 0) {
				for (int j = i + 1; j < Set_size; j++) {
					a2 = solutionSet.get(j);
					int num1 = 0, num3 = 0;
					for (int k = 0; k < num_obj; k++) {
						if (a1.getObjective(k) > a2.getObjective(k)) {
							num1++;
						} else if (a1.getObjective(k) < a2.getObjective(k)) {
							num3++;
						}
						//else
						//	num3++;
					}
					if (num1 == 0) {
						dominateMe[j] = 1;
					} else if (num3 == 0) {
						dominateMe[i] = 1;
						break;
					}
				}
				if (dominateMe[i] == 0) {
					ranking_.add(solutionSet.get(i));
				}
			}
		}
		return ranking_;

	} // Ranking

	/**
	 * Runs of the SMPSO algorithm.
	 *
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 * solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initParams();

		// ->Step 1 (and 3) Create the initial population and evaluate
		for (int i = 0; i < all_swarmSize_; i++) {
			Solution particle = new Solution(problem_);
			problem_.evaluate(particle);
			problem_.evaluateConstraints(particle);
			particles_.add(particle);
		}

		// -> Step2. Randomly Initialize the speed_ of each particle 
		for (int i = 0; i < all_swarmSize_; i++) {
			for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
				double r1 = PseudoRandom.randDouble(0, 1);
				speed_[i][j] = r1 * deltaMax_[j];
				r1 = PseudoRandom.randDouble(0, 1);
				if (r1 < 0.5) {
					speed_[i][j] = -speed_[i][j];
				}
				//speed_[i][j]=0.0;
			}
		}

		//// -> Step 3. update the pbest of each particle
		for (int i = 0; i < particles_.size(); i++) {
			Solution particle = new Solution(particles_.get(i));
			pbest_[i] = particle;
		}

		//// -> Step 4. update the gbest of each swarm
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			int maxnum = swarmSize_ * (i + 1);
			int gbest_index = -1;
			double fitness = Double.MAX_VALUE;
			for (int j = i * swarmSize_; j < maxnum; j++) {
				if (fitness > particles_.get(j).getObjective(i)) {
					gbest_index = j;
					fitness = particles_.get(j).getObjective(i);
				}
			}
			Solution particle = new Solution(particles_.get(gbest_index));
			gbest_[i] = particle;
		}

		// Step 5 update the archive 
		for (int i = 0; i < particles_.size(); i++) {
			Solution particle = new Solution(particles_.get(i));
			archive_.add(particle);
		}
		archive_ = Getnondominated(archive_);


		// -> Step 7. Iterations ..
		Random r = new Random();
		function_ev = 0;
		while (function_ev < max_function_ev) {
			try {
				// Compute the speed_
				computeSpeed(function_ev, max_function_ev);
			} catch (IOException ex) {
				Logger.getLogger(SMPSO.class.getName()).log(Level.SEVERE, null,
						ex);
			}

			// Compute the new positions for the particles_
			computeNewPositions();

			// Evaluate the new particles_ in new positions
			for (int i = 0; i < particles_.size(); i++) {
				Solution particle = particles_.get(i);
				problem_.evaluate(particle);
				problem_.evaluateConstraints(particle);
			}

			// Actualize the  of this particle
			for (int i = 0; i < particles_.size(); i++) {
				int j = i / swarmSize_;
				if (particles_.get(i).getObjective(j) < pbest_[i].getObjective(j)) {
					Solution particle = new Solution(particles_.get(i));
					pbest_[i] = particle;
					if (pbest_[i].getObjective(j) < gbest_[j].getObjective(j)) {
						particle = new Solution(pbest_[i]);
						gbest_[j] = particle;
					}

				}
			}

			int sizeofarchive = archive_.size();
			//perform ELS on the original particle in the archive
			for (int i = 0; i < sizeofarchive; i++) {
				Solution one_particle = new Solution(archive_.get(i));
				XReal particle = new XReal(one_particle);
				int seleted_ind;
				double value;
				seleted_ind = PseudoRandom.randInt(0, problem_.getNumberOfVariables() - 1);
				value = 0.25 * r.nextGaussian();
				particle.setValue(seleted_ind, particle.getValue(seleted_ind) +
						value * (problem_.getUpperLimit(seleted_ind) - problem_.getLowerLimit(seleted_ind)));
				if (particle.getValue(seleted_ind) < problem_.getLowerLimit(seleted_ind)) {
					particle.setValue(seleted_ind, problem_.getLowerLimit(seleted_ind));
				}
				if (particle.getValue(seleted_ind) > problem_.getUpperLimit(seleted_ind)) {
					particle.setValue(seleted_ind, problem_.getUpperLimit(seleted_ind));
				}
				problem_.evaluate(one_particle);
				problem_.evaluateConstraints(one_particle);
				archive_.add(one_particle);
			}

			// add the new pbest into archive
			for (int i = 0; i < particles_.size(); i++) {
				Solution particle = new Solution(pbest_[i]);
				archive_.add(particle);
			}

			archive_ = Getnondominated(archive_);//get the nondominated solutions from archive

			// Assign crowding distance to the leaders_
			if (archive_.size() > archiveSize_) {
				distance_.crowdingDistanceAssignment(archive_,
						problem_.getNumberOfObjectives());
				archive_.sort(new CrowdingComparator());
				while (archive_.size() > archiveSize_) {
					archive_.remove(archive_.size() - 1);
				}
			}
			function_ev += particles_.size();
			function_ev += sizeofarchive;
			//iteration_++;
			/*
			 * if ((iteration_ % 1) == 0) {
			 * leaders_.printObjectivesOfValidSolutionsToFile("FUNV"+iteration_)
			 * ; leaders_.printObjectivesToFile("FUN"+iteration_) ;
			 * leaders_.printVariablesToFile("VAR"+iteration_) ; }
			 */
		}
		// leaders_.printObjectivesOfValidSolutionsToFile("FUNV.SMPSO") ;
		archive_.printFeasibleFUN("FUN_SMPSO");

		return this.archive_;
	} // execute

} // SMPSO