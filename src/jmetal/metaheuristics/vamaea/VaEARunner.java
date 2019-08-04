/**
 * RaPSO.java
 *
 * @author Xin.Hu
 */
package jmetal.metaheuristics.vamaea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgHvPlatEMO;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.plot.LineBeyend4d;
import jmetal.util.plot.Scatter2d;
import jmetal.util.plot.Scatter3d;
import org.jfree.ui.RefineryUtilities;

import javax.swing.*;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;


public class VaEARunner {

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException {
		// the numbers of objectives
		int m = 6;
		final int low = 6;
		Logger logger = Configuration.logger_;
		FileHandler fileHandler = new FileHandler("r2pso.log");
		logger.addHandler(fileHandler);
		for (int fun = low; fun <= low; fun++) {
			Problem problem = null;
			Algorithm algorithm;
			QualityIndicator indicators;
			indicators = null;
			boolean wfgIs2d = false;
			if (args.length == 1) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem
				problem = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getindicator();
			}
			// init parameter of algorithm
			algorithm = new VaEA(problem);
			coffientSetting(algorithm, problem, fun);
			SolutionSet population;
			long initTime = System.currentTimeMillis();
			population = algorithm.execute();
			long endTime = System.currentTimeMillis() - initTime;
			plot(problem, population, indicators, true);
			logger.info("Total run time is" + endTime + "ms");
			wfgHvPlatEMO wfgHvPlatEMO = new wfgHvPlatEMO(population.writeObjectivesToMatrix(), indicators.getTrueParetoFront());
			double hv = wfgHvPlatEMO.calculatewfghv();
			assert indicators != null;
			logger.info(problem.getName()
					+ "\nHyperVolume: " + hv
					+ "\nEPSILON    : " + indicators.getEpsilon(population)
					+ "\nGD         : " + indicators.getGD(population)
					+ "\nIGD        : " + indicators.getCEC_IGD(population)
					+ "\nSpread     : " + indicators.getGeneralizedSpread(population)
					+ "\nSpace        : " + indicators.getSpace(population)
					+ "\nNumberOfPF        : " + population.size()
					+ "\nPD                : " + indicators.getPD(population)
			);
		}
	}

	private static void coffientSetting(Algorithm algorithm, Problem problem, int fun) throws JMException {
		Operator clone = null;
		// Crossover operator
		Operator crossover;
		// Mutation operator
		Operator mutation;
		if (problem.getNumberOfObjectives() == 2) {
			if (fun < 6) {
				algorithm.setInputParameter("maxIterations", 250);
			} else if (fun < 22) {
				algorithm.setInputParameter("maxIterations", 500);
			} else {
				algorithm.setInputParameter("maxIterations", 1000);
			}
			algorithm.setInputParameter("swarmSize", 100);
			HashMap<String, Integer> parameters = new HashMap<String, Integer>();
			parameters.put("clonesize", 100);
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
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
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
		} else if (problem.getNumberOfObjectives() == 5) {
			algorithm.setInputParameter("maxIterations", 500);
			algorithm.setInputParameter("swarmSize", 126);
			// Clone operator
			HashMap<String, Integer> parameters = new HashMap<String, Integer>();
			parameters.put("clonesize", 210);
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
		} else if (problem.getNumberOfObjectives() == 6) {
			algorithm.setInputParameter("maxIterations", 500);
			algorithm.setInputParameter("swarmSize", 132);
			// Clone operator
			HashMap<String, Integer> parameters = new HashMap<String, Integer>();
			parameters.put("clonesize", 132);
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
		} else if (problem.getNumberOfObjectives() == 8) {
			algorithm.setInputParameter("maxIterations", 1000);
			algorithm.setInputParameter("swarmSize", 156);
			// Clone operator
			HashMap<String, Integer> parameters = new HashMap<String, Integer>();
			parameters.put("clonesize", 156);
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
		} else if (problem.getNumberOfObjectives() == 10) {
			algorithm.setInputParameter("maxIterations", 500);
			algorithm.setInputParameter("swarmSize", 275);
			// Clone operator
			HashMap<String, Integer> parameters = new HashMap<String, Integer>();
			parameters.put("clonesize", 275);
			clone = CloneFactory.getClone("ShiftedDistanceClone", parameters);
		}
		HashMap<String, Double> parameters = new HashMap<String, Double>();
//		parameters.put("CR", 0.2);
//		parameters.put("F", 0.5);
//		crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

		parameters.put("probability", 1.0);
		parameters.put("distributionIndex", 30.0);
		crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
		parameters = new HashMap<>();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		parameters.put("distributionIndex", 20.0);
		mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
		parameters = null;
		Operator selection;
		selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);//

		algorithm.setInputParameter("normalize", true);
		// Add the operators to the algorithm
		algorithm.addOperator("clone", clone);
		algorithm.addOperator("crossover", crossover);
		algorithm.addOperator("mutation", mutation);
		algorithm.addOperator("selection", selection);

	}

	public static void plot(Problem problem, SolutionSet population, QualityIndicator indicators, boolean truePF) {
		if (2 == problem.getNumberOfObjectives()) {
			final Scatter2d demo = new Scatter2d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, truePF);
			demo.pack();
			RefineryUtilities.centerFrameOnScreen(demo);
			demo.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
			demo.setSize(1000, 720);
			demo.setVisible(true);
		} else if (3 == problem.getNumberOfObjectives()) {
			new Scatter3d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, truePF).plot();
		} else {
			final LineBeyend4d demo = new LineBeyend4d("Dimension", "Fitness", problem.getName(), population.writeObjectivesToMatrix(), indicators);
			demo.pack();
			RefineryUtilities.centerFrameOnScreen(demo);
			demo.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
			demo.setSize(1280, 720);
			demo.setVisible(true);
		}
	}
}