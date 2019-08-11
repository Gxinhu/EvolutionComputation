package jmetal.metaheuristics.nsgaIII;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.WFG.WFG5;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgHvPlatEMO;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.HashMap;

public class NSGAIII_main {
	public static void main(String args[]) throws JMException, ClassNotFoundException, IOException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; //Selection operator

		HashMap parameters; // Operator parameters

		problem = new WFG5("Real", 14, 20, 8);
		algorithm = new NSGAIII(problem);


		algorithm.setInputParameter("normalize", true);


		algorithm.setInputParameter("div1", 3);
		algorithm.setInputParameter("div2", 2);


		algorithm.setInputParameter("maxGenerations", 500);

		// Mutation and Crossover for Real codification
		parameters = new HashMap();
		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 30.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
				parameters);

		parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation",
				parameters);

		parameters = null;
		selection = SelectionFactory.getSelectionOperator("RandomSelection",
				parameters);

		// Add the operators to the algorithm
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);

		SolutionSet population = algorithm.execute();
		wfgHvPlatEMO wfgHvPlatEMO = new wfgHvPlatEMO(population.writeObjectivesToMatrix(), problem.getName());
		double hv1 = wfgHvPlatEMO.calculatewfghv();
		System.out.println(hv1);
	}
}
