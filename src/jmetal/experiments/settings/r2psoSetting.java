//  NSGAIISettings.java

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.r2pso.VagPso;
import jmetal.operators.clone.Clone;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.Crossover;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.Mutation;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.util.JMException;

import java.util.HashMap;

/**
 * Settings class of algorithm PaRPEA (real encoding)
 */
public class r2psoSetting extends Settings {
	public int populationSize_;
	public int maxEvaluations_;
	public int maxGenerations_;
	public double mutationProbability_;
	public double crossoverProbability_;
	public double mutationDistributionIndex_;
	public double crossoverDistributionIndex_;

	// For Permutation variable
	public r2psoSetting(String problem, Object[] params) {
		super(problem);
		try {
			problem_ = (new ProblemFactory()).getProblem(problemName_, params);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Default experiments.settings
		if (problem_.getNumberOfObjectives() == 2) {

			if (problem_.getName().contains("UF")
					|| problem_.getName().contains("LZ09")) {
				populationSize_ = 100;
				maxGenerations_ = 3000;
			} else if (problem_.getName().contains("ZHX")) {
				populationSize_ = 80;
				maxGenerations_ = 375;
			} else {
				populationSize_ = 80;
				maxGenerations_ = 500;
			}

		}
		if (problem_.getNumberOfObjectives() == 3) {
			populationSize_ = 105;

			maxGenerations_ = 500;

		} else if (problem_.getNumberOfObjectives() == 5) {

			populationSize_ = 126;
			maxGenerations_ = 200;
		} else if (problem_.getNumberOfObjectives() == 6) {

			populationSize_ = 132;
			maxGenerations_ = 500;
		} else if (problem_.getNumberOfObjectives() == 8) {
			populationSize_ = 156;
			maxGenerations_ = 400;

		} else if (problem_.getNumberOfObjectives() == 9) {

			populationSize_ = 210;
			maxGenerations_ = 2000;

		} else if (problem_.getNumberOfObjectives() == 10) {

			populationSize_ = 275;
			maxGenerations_ = 300;

		} else if (problem_.getNumberOfObjectives() == 15) {
			populationSize_ = 135;
			maxGenerations_ = 500;
		} else if (problem_.getNumberOfObjectives() == 20) {
			populationSize_ = 230;
			maxGenerations_ = 250;
		}
//		if (problem_.getNumberOfObjectives() <= 10) {
//			maxGenerations_ = 100000 / populationSize_;
//		} else {
//			maxGenerations_ = 10000 * problem_.getNumberOfObjectives() / populationSize_;
//		}
		crossoverProbability_ = 1.0;
		mutationProbability_ = 1.0 / problem_.getNumberOfVariables();
		crossoverDistributionIndex_ = 30.0;
		mutationDistributionIndex_ = 20.0;
		// NSGAIISettings
	}

	/**
	 * Configure AdMOEA with user-defined parameter experiments.settings
	 *
	 * @return A AdMOEA algorithm object
	 * @throws JMException
	 */
	@Override
	public Algorithm configure() throws JMException {
		Algorithm algorithm;
		Crossover crossover;
		Mutation mutation;
		Clone clone;

		HashMap parameters; // Operator parameters

		// Creating the algorithm.

		algorithm = new VagPso(problem_);
		// Algorithm parameters
		algorithm.setInputParameter("maxIterations", maxGenerations_);
		algorithm.setInputParameter("swarmSize", populationSize_);

		/**
		 * Mutation and Crossover for Real codification
		 */
		parameters = new HashMap();
		parameters.put("probability", crossoverProbability_);
		parameters.put("distributionIndex", crossoverDistributionIndex_);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
//		parameters = new HashMap();
//		parameters.put("CR", 0.2);
//		parameters.put("F", 0.5);
//		crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

		parameters.put("probability", mutationProbability_);
		parameters.put("distributionIndex", mutationDistributionIndex_);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
		parameters.put("clonesize", populationSize_);
		clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
//		// mutation = MutationFactory.getMutationOperator("GaussMutation",
//		// parameters);
		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("clone", clone);
		return algorithm;
	} // configure
} // NSGAIISettings
