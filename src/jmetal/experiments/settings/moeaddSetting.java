//  NSGAIISettings.java

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.moeadd.MOEADD;
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
public class moeaddSetting extends Settings {
	public int populationSize_;
	public int maxEvaluations_;
	public int maxGenerations;
	public double mutationProbability_;
	public double crossoverProbability_;
	public double mutationDistributionIndex_;
	public double crossoverDistributionIndex_;

	// For Permutation variable
	public moeaddSetting(String problem, int populationSize_, int maxGenerations, Object[] params) {
		super(problem);
		try {
			problem_ = (new ProblemFactory()).getProblem(problemName_, params);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Default experiments.settings
		this.populationSize_ = populationSize_;
		this.maxGenerations = maxGenerations;
		maxEvaluations_ = this.maxGenerations * this.populationSize_;
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

		HashMap parameters; // Operator parameters

		// Creating the algorithm.
		algorithm = new MOEADD(problem_);

		// Algorithm parameters
		algorithm.setInputParameter("maxEvaluations", maxEvaluations_);
		algorithm.setInputParameter("populationSize", populationSize_);
		algorithm.setInputParameter("dataDirectory", "weight");

		/**
		 * Mutation and Crossover for Real codification
		 */
		parameters = new HashMap();
		parameters.put("probability", crossoverProbability_);
		parameters.put("distributionIndex", crossoverDistributionIndex_);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
		parameters = new HashMap();
		parameters.put("probability", mutationProbability_);
		parameters.put("distributionIndex", mutationDistributionIndex_);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		return algorithm;
	} // configure
} // NSGAIISettings
