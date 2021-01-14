//  NSGAIISettings.java

package jmetal.experiments.settings;

import jmetal.core.Algorithm;
import jmetal.experiments.Settings;
import jmetal.metaheuristics.Vapso.VaPSO;
import jmetal.metaheuristics.Vapso.VagPsoconstant;
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
public class VaPSOSetting extends Settings {
	public int populationSize_;
	public int maxEvaluations_;
	public int maxGenerations_;
	public double mutationProbability_;
	public double crossoverProbability_;
	public double mutationDistributionIndex_;
	public double crossoverDistributionIndex_;

	public VaPSOSetting(String problem, int populationSize_, int maxGenerations, Object[] params) {
		super(problem);
		try {
			problem_ = (new ProblemFactory()).getProblem(problemName_, params);
		} catch (JMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Default experiments.settings
		this.populationSize_ = populationSize_;
		this.maxGenerations_ = maxGenerations;
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

		algorithm = new VaPSO(problem_);
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
