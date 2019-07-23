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
	public moeaddSetting(String problem, Object[] params) {
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
				maxGenerations = 3000;
			} else if (problem_.getName().contains("ZHX")) {
				populationSize_ = 80;
				maxGenerations = 375;
			} else {
				populationSize_ = 80;
				maxGenerations = 500;
			}

		}
		if (problem_.getNumberOfObjectives() == 3) {
			populationSize_ = 105;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 500;
			} else {
				maxGenerations = 500;
			}
//			maxGenerations = 1000;
			if (problem_.getName().contains("UF")
					|| problem_.getName().contains("LZ09")) {
				populationSize_ = 300;
				maxGenerations = 1000;
			}

		} else if (problem_.getNumberOfObjectives() == 5) {

			populationSize_ = 212;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 500;
			} else {
				maxGenerations = 500;
			}
		} else if (problem_.getNumberOfObjectives() == 6) {

			populationSize_ = 132;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 500;
			} else {
				maxGenerations = 500;
			}
		} else if (problem_.getNumberOfObjectives() == 8) {
			populationSize_ = 156;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 500;
			} else {
				maxGenerations = 500;
			}

		} else if (problem_.getNumberOfObjectives() == 9) {

			populationSize_ = 210;
			maxGenerations = 2000;

		} else if (problem_.getNumberOfObjectives() == 10) {

			populationSize_ = 275;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 500;
			} else {
				maxGenerations = 500;
			}

		}
		maxEvaluations_ = maxGenerations * populationSize_;

		crossoverProbability_ = 1.0;
		mutationProbability_ = 1.0 / problem_.getNumberOfVariables();
		crossoverDistributionIndex_ = 20.0;
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
