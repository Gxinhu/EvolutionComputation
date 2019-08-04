/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */
package jmetal.metaheuristics.agmopso;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;


public class AgMOPSO_main {
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		long startime = System.currentTimeMillis();
		int m = 6;
		int low = 1;
		int high = 2;
		System.out.println("IGDmean" + "\t" + "std" + "\t" + "GDmean" + "\t" + "std" + "\t" + "GSpreadmean" + "\t" + "std" + "\t" + "problemName");
		for (int fun = low; fun <= high; fun++) {
			int runtimes = 1;
			double[] IGDarray = new double[runtimes];
			double[] GDArray = new double[runtimes];
			double[] GeneraltionalSpread = new double[runtimes];
			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator clone = null; // Crossover operator
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			QualityIndicator indicators; // Object to get quality indicators
			indicators = null;
			boolean wfgis2d = true;
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
				problem = new cricleselectproblem(problem, indicators, fun, m, wfgis2d).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgis2d).getindicator();
			}

//                algorithm = new AgmopsoAdaptiveWeighter(problem, indicators, i);
//                algorithm=new AgR2ADW(problem,indicators,i);
//				algorithm = new AgmopsoDE(problem, indicators, i, false);
//				algorithm = new AgMOPSOwithR2newVersion(problem);
			int k = 0;
			algorithm = new AgMOPSO(problem, indicators, k, false);
			if (problem.getNumberOfObjectives() == 2) {
				if (fun < 6) {
					algorithm.setInputParameter("maxIterations", 250);
				} else if (fun < 22) {
					algorithm.setInputParameter("maxIterations", 500);
				} else {
					algorithm.setInputParameter("maxIterations", 3000);
				}
				algorithm.setInputParameter("swarmSize", 100);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 100);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 3) {
				if (fun < 22) {
					algorithm.setInputParameter("maxIterations", 500);
				} else {
					algorithm.setInputParameter("maxIterations", 3000);
				}
				algorithm.setInputParameter("swarmSize", 105);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 105);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 5) {
				algorithm.setInputParameter("maxIterations", 500);
				algorithm.setInputParameter("swarmSize", 210);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 210);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 6) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 132);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 132);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 8) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 156);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 156);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 10) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 275);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 275);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			}
			HashMap<String, Double> parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0);
			parameters.put("distributionIndex", 20.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
//				parameters.put("CR", 1.0);
//				parameters.put("F", 0.5);
//				crossover = CrossoverFactory.getCrossoverOperator(
//						"DifferentialEvolutionCrossover", parameters);
			parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			// Add the operators to the algorithm
			SolutionSet population = null;
			algorithm.addOperator("clone", clone);
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			for (int i = 0; i < runtimes; i++) {
				population = algorithm.execute();
				IGDarray[i] = indicators.getCEC_IGD(population);
				GeneraltionalSpread[i] = indicators.getGeneralizedSpread(population);
				GDArray[i] = indicators.getGD(population);
			}
			population.printObjectivesToFile("./" + fun + ".txt");
			TestStatistics sta = null;
			Arrays.sort(IGDarray);
			sta = new TestStatistics(IGDarray);
			System.out.print(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t");
			Arrays.sort(GDArray);
			sta = new TestStatistics(GDArray);
			System.out.print(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t");
			Arrays.sort(GeneraltionalSpread);
			sta = new TestStatistics(GeneraltionalSpread);
			System.out.println(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getNumberOfObjectives() + problem.getName());

		} //runtimes
	} // main
}// AgMOPSO_main