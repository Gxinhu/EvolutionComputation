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
	public NSGAIIISetting(String problem, int populationSize, int generations, Object[] params) {
		super(problem);
		try {
			problem_ = (new ProblemFactory()).getProblem(problemName_, params);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Default experiments.settings
		this.populationSize = populationSize;
		this.maxGenerations = generations;
		maxEvaluations = this.maxGenerations * this.populationSize;
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
