package jmetal.metaheuristics.ragpso;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.agmopso.TestStatistics;
import jmetal.metaheuristics.selectProblemRvea;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculateRvea;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

public class ragmopsoDE_SBXmain {
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException {
		// the number of objectives
		int m = 3;
		final int low = 6;
		final int high = 8;
		final int runtime = 30;
		double[] hv = new double[runtime];
		for (int fun = low; fun <= high; fun++) {
			// The problem to solve
			Problem problem = null;
			// The algorithm to use
			Algorithm algorithm;
			// Crossover operator
			Operator crossoverDe;
			Operator crossoverSbx;
			// Mutation operator
			Operator mutation;
			// Selection operator
			Operator selection;
			// Object to get quality indicators
			QualityIndicator indicators = null;
			//choose the problem
			if (args.length == 1) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = {"Real"};
				indicators = new QualityIndicator(problem, args[1]);
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else { // Default problem
				problem = new selectProblemRvea(problem, indicators, fun, m).getProblem();
				indicators = new selectProblemRvea(problem, indicators, fun, m).getindicator();
			}
			// init parameter of algorithm
			int k = 0;
			algorithm = new ragmopsoVersion4DE_SBX(problem, indicators, k);

			if (fun == 6 | fun == 8) {
				algorithm.setInputParameter("maxIterations", 1000);
			} else if (fun <= 12) {
				algorithm.setInputParameter("maxIterations", 500);
			} else if (fun > 12 & fun < 22) {
				algorithm.setInputParameter("maxIterations", 1000);
			}
			if (problem.getNumberOfObjectives() == 2) {
				algorithm.setInputParameter("swarmSize", 100);
			} else if (problem.getNumberOfObjectives() == 3) {
				algorithm.setInputParameter("swarmSize", 105);
			} else if (problem.getNumberOfObjectives() == 5) {
				algorithm.setInputParameter("swarmSize", 210);
			} else if (problem.getNumberOfObjectives() == 6) {
				algorithm.setInputParameter("swarmSize", 132);
			} else if (problem.getNumberOfObjectives() == 8) {
				algorithm.setInputParameter("swarmSize", 156);
			} else if (problem.getNumberOfObjectives() == 10) {
				algorithm.setInputParameter("swarmSize", 275);
			}
			HashMap<String, Double> parameters = new HashMap<>();
			parameters.put("probability", 1.0);
			parameters.put("distributionIndex", 30.0);
			crossoverSbx = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
			parameters.put("CR", 0.1);
			parameters.put("F", 0.5);
			crossoverDe = CrossoverFactory.getCrossoverOperator(
					"DifferentialEvolutionCrossover", parameters);
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
			selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection",
					null);
			algorithm.addOperator("selection", selection);
			// Add the operators to the algorithm
			SolutionSet population;
			algorithm.addOperator("crossoverDe", crossoverDe);
			algorithm.addOperator("crossoverSbx", crossoverSbx);
			algorithm.addOperator("mutation", mutation);
			for (int i = 0; i < runtime; i++) {
				population = algorithm.execute();
				wfghvCalculateRvea wfg = new wfghvCalculateRvea(population, fun);
				hv[i] = wfg.calculatewfghv();
			}
			Arrays.sort(hv);
			TestStatistics sta = new TestStatistics(hv);
			System.out.println(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getNumberOfObjectives() + problem.getName());

		}

	} // main
}// AgMOPSO_main