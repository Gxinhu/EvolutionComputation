/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */
package jmetal.metaheuristics.Vapso;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.agmopso.TestStatistics;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgHvPlatEMO;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;


public class r2psoUseShiftedDistanceMain {
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		long startime = System.currentTimeMillis();
		int m = 3;
		int low = 10;
		int high = 10;
		System.out.println("IGDmean" + "\t" + "std" + "\t" + "GDmean" + "\t" + "std" + "\t" + "GSpreadmean" + "\t" + "std" + "\t" + "problemName");
		for (int fun = low; fun <= high; fun++) {
			int runTimes = 30;
			double[] iGD = new double[runTimes];
			double[] gD = new double[runTimes];
			double[] generaltionalSpread = new double[runTimes];
			// The problem to solve
			Problem problem = null;
			/* The algorithm to use */
			Algorithm algorithm;
			// Crossover operator
			Operator clone = null;
			// Crossover operator
			Operator crossover;
			// Mutation operator
			Operator mutation;
			// Object to get quality indicators
			QualityIndicator indicators;
			indicators = null;
			boolean wfgis2d = false;
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


			algorithm = new VagPsoAdaptivelbest(problem);
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
				algorithm.setInputParameter("maxIterations", 100000 / 252);
				algorithm.setInputParameter("swarmSize", 252);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 252);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 8) {
				algorithm.setInputParameter("maxIterations", 100000 / 330);
				algorithm.setInputParameter("swarmSize", 330);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 330);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 10) {
				algorithm.setInputParameter("maxIterations", 300);
				algorithm.setInputParameter("swarmSize", 275);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 275);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 15) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 135);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 135);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			}
			HashMap<String, Double> parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0);
			parameters.put("distributionIndex", 30.0);
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
			double startTime = System.currentTimeMillis();
			double[] hv = new double[runTimes];
			for (int i = 0; i < runTimes; i++) {
				population = algorithm.execute();
				wfgHvPlatEMO wfg = new wfgHvPlatEMO(population.writeObjectivesToMatrix(), problem.getName());
				iGD[i] = indicators.getCEC_IGD(population);
				generaltionalSpread[i] = indicators.getGeneralizedSpread(population);
				gD[i] = indicators.getGD(population);
				hv[i] = wfg.calculatewfghv();
			}
			double endTime = System.currentTimeMillis();
			TestStatistics sta;
			Arrays.sort(iGD);
			sta = new TestStatistics(iGD);
			System.out.println("\t" + sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getName());
//			Arrays.sort(gD);
//			sta = new TestStatistics(gD);
//			System.out.print(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t");
//			Arrays.sort(generaltionalSpread);
//			sta = new TestStatistics(generaltionalSpread);
//			System.out.println(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getNumberOfObjectives() + problem.getName());
			Arrays.sort(hv);
			sta = new TestStatistics(hv);
			System.out.println("\t" + sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getName());
		} //runtimes
	} // main
}// AgMOPSO_main