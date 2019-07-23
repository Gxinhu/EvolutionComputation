//  NSGAIISettings.java

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.NSGAIII.NSGAIII_SBX;
import jmetal.operators.crossover.Crossover;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.Mutation;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.Selection;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.util.JMException;

import java.util.HashMap;

/**
 * Settings class of algorithm PaRPEA (real encoding)
 */
public class NSGAIIISetting extends Settings {
	public int populationSize;
	public int maxEvaluations;
	public int maxGenerations;
	public double mutationProbability;
	public double crossoverProbability;
	public double mutationDistributionIndex;
	public double crossoverDistributionIndex;
	public int div1;
	private int div2;

	// For Permutation variable
	public NSGAIIISetting(String problem, Object[] params) {
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
				populationSize = 100;
				maxGenerations = 3000;
			} else if (problem_.getName().contains("ZHX")) {
				populationSize = 80;
				maxGenerations = 375;
			} else {
				populationSize = 80;
				maxGenerations = 500;
			}

		}
		if (problem_.getNumberOfObjectives() == 3) {
			populationSize = 105;
			div1 = 13;
			div2 = 0;
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
				populationSize = 300;
				maxGenerations = 1000;
			}

		} else if (problem_.getNumberOfObjectives() == 5) {

			populationSize = 212;
			div1 = 6;
			div2 = 0;
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
				maxGenerations = 1000;
			}

		} else if (problem_.getNumberOfObjectives() == 6) {
			populationSize = 132;
			div1 = 4;
			div2 = 1;
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
			populationSize = 156;
			div1 = 3;
			div2 = 2;
			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 400;
			} else {
				maxGenerations = 500;
			}

		} else if (problem_.getNumberOfObjectives() == 9) {

			populationSize = 210;
			maxGenerations = 2000;

		} else if (problem_.getNumberOfObjectives() == 10) {

			populationSize = 276;
			div1 = 3;
			div2 = 2;
			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ5")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ6")) {
				maxGenerations = 300;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ7")) {
				maxGenerations = 300;
			} else {
				maxGenerations = 500;
			}

		} else if (problem_.getNumberOfObjectives() == 15) {

			populationSize = 136;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 3000;
			} else {
				maxGenerations = 3000;
			}
		} else if (problem_.getNumberOfObjectives() == 20) {

			populationSize = 200;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 3000;
			} else {
				maxGenerations = 3000;
			}
		} else if (problem_.getNumberOfObjectives() == 25) {

			populationSize = 300;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations = 3000;
			} else {
				maxGenerations = 3000;
			}
		}

		maxEvaluations = maxGenerations * populationSize;
		crossoverProbability = 1.0;
		mutationProbability = 1.0 / problem_.getNumberOfVariables();
		crossoverDistributionIndex = 30.0;
		mutationDistributionIndex = 20.0;
	} // NSGAIISettings

	/**
	 * Configure AdMOEA with user-defined parameter experiments.settings
	 *
	 * @return A AdMOEA algorithm object
	 * @throws JMException
	 */
	@Override
	public Algorithm configure() throws JMException {
		Algorithm algorithm;
		Selection selection;
		Crossover crossover;
		Mutation mutation;

		HashMap parameters; // Operator parameters

		// Creating the algorithm.
		algorithm = new NSGAIII_SBX(problem_);

		// Algorithm parameters
		algorithm.setInputParameter("maxEvaluations", maxEvaluations);
		algorithm.setInputParameter("swarmSize", populationSize);
		algorithm.setInputParameter("normalize", true);

		/**
		 * Mutation and Crossover for Real codification
		 */
		parameters = new HashMap();
		parameters.put("probability", crossoverProbability);
		parameters.put("distributionIndex", crossoverDistributionIndex);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);


		algorithm.setInputParameter("div1", div1);
		algorithm.setInputParameter("div2", div2);
		parameters.put("probability", mutationProbability);
		parameters.put("distributionIndex", mutationDistributionIndex);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
		parameters = null;
		selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);
		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		return algorithm;
	} // configure
} // NSGAIISettings
