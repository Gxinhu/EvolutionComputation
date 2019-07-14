//  NSGAIISettings.java

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.nsgaII.NSGAII;
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
public class NSGAIISettings extends Settings {
	public int populationSize_;
	public int maxEvaluations_;
	public int maxGenerations_;
	public double mutationProbability_;
	public double crossoverProbability_;
	public double mutationDistributionIndex_;
	public double crossoverDistributionIndex_;

	// For Permutation variable
	public NSGAIISettings(String problem, Object[] params) {
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
			populationSize_ = 92;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 400;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 250;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 600;
			} else {
				maxGenerations_ = 1000;
			}
//			maxGenerations_ = 1000;
			if (problem_.getName().contains("UF")
					|| problem_.getName().contains("LZ09")) {
				populationSize_ = 300;
				maxGenerations_ = 1000;
			}

		} else if (problem_.getNumberOfObjectives() == 5) {

			populationSize_ = 212;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 600;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 350;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 1000;
			} else {
				maxGenerations_ = 1250;
			}

		} else if (problem_.getNumberOfObjectives() == 8) {
			populationSize_ = 156;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 750;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 1250;
			} else {
				maxGenerations_ = 1500;
			}

		} else if (problem_.getNumberOfObjectives() == 9) {

			populationSize_ = 210;
			maxGenerations_ = 2000;

		} else if (problem_.getNumberOfObjectives() == 10) {

			populationSize_ = 276;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 750;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 2000;
			} else {
				maxGenerations_ = 2000;
			}

		} else if (problem_.getNumberOfObjectives() == 15) {

			populationSize_ = 136;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 3000;
			} else {
				maxGenerations_ = 3000;
			}
		} else if (problem_.getNumberOfObjectives() == 20) {

			populationSize_ = 200;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 3000;
			} else {
				maxGenerations_ = 3000;
			}
		} else if (problem_.getNumberOfObjectives() == 25) {

			populationSize_ = 300;

			if (problem_.getName().equalsIgnoreCase("DTLZ1")) {
				maxGenerations_ = 1500;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ2")) {
				maxGenerations_ = 1000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ3")) {
				maxGenerations_ = 2000;
			} else if (problem_.getName().equalsIgnoreCase("DTLZ4")) {
				maxGenerations_ = 3000;
			} else {
				maxGenerations_ = 3000;
			}
		}

//		maxGenerations_ = maxGenerations_/2;
		maxEvaluations_ = maxGenerations_ * populationSize_;

//		maxEvaluations_ = 25000; // Use the same evaluations 
		crossoverProbability_ = 1.0;
		mutationProbability_ = 1.0 / problem_.getNumberOfVariables();
		crossoverDistributionIndex_ = 20.0;
		mutationDistributionIndex_ = 20.0;
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
		algorithm = new NSGAII(problem_, false, 0); //

		// Algorithm parameters
		algorithm.setInputParameter("maxIterations", maxGenerations_);
		algorithm.setInputParameter("maxEvaluations", maxEvaluations_);
		algorithm.setInputParameter("swarmSize", populationSize_);

		/**
		 * Mutation and Crossover for Real codification
		 */
		parameters = new HashMap();
		parameters.put("probability", crossoverProbability_);
		parameters.put("distributionIndex", crossoverDistributionIndex_);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

		// Crossover operator
//		 parameters = new HashMap() ;
//		 parameters.put("CR", 1.0) ;
//		 parameters.put("F", 0.5) ;
//		 parameters.put("DE_VARIANT", "rand/1/exp") ;//rand/1/binrand/1/exp
//		 crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",
//		 parameters);

		parameters = new HashMap();
		parameters.put("probability", mutationProbability_);
		parameters.put("distributionIndex", mutationDistributionIndex_);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
//		// mutation = MutationFactory.getMutationOperator("GaussMutation",
//		// parameters);
		parameters = null;
		selection = SelectionFactory.getSelectionOperator("BinaryTournament",
				parameters);
		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);
		return algorithm;
	} // configure
} // NSGAIISettings
