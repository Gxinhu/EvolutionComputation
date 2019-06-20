package jmetal.metaheuristics.ragpso;
/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.selectProblemRvea;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculateRvea;
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


public class ragmopsorunnerDE {
	public static Logger logger_; // Logger object


	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException {
		// the number of objectives
		int m = 6;
		logger_ = Configuration.logger_;
		// FileHandler object
		FileHandler fileHandler_ = new FileHandler("ragmopso.log");
		logger_.addHandler(fileHandler_);
		final int low = 8;
		for (int fun = low; fun <= low; fun++) {
			// The problem to solve
			Problem problem = null;
			// The algorithm to use
			Algorithm algorithm;
			// Crossover operator
			Operator crossover;
			// Mutation operator
			Operator mutation;
			// Selection operator
			Operator selection;
			QualityIndicator indicators; // Object to get quality indicators
			indicators = null;
			//choose the problem
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
				problem = new selectProblemRvea(problem, indicators, fun, m).getProblem();
				indicators = new selectProblemRvea(problem, indicators, fun, m).getindicator();
			}
			// init parameter of algorithm
			int i = 0;
			algorithm = new ragmopsoversion4DE(problem, indicators, i);

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

			// init operator
			HashMap<String, Double> parameters = new HashMap<>();
			parameters.put("CR", 0.1);
			parameters.put("F", 0.5);
			crossover = CrossoverFactory.getCrossoverOperator(
					"DifferentialEvolutionCrossover", parameters);
			parameters = new HashMap<>();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			parameters = null;
			selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection",
					parameters);
			algorithm.addOperator("selection", selection);

			// Add the operators to the algorithm
			SolutionSet population;
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			long initTime = System.currentTimeMillis();
			population = algorithm.execute();
			long endTime = System.currentTimeMillis() - initTime;

			// plot
			if (2 == problem.getNumberOfObjectives()) {
				final Scatter2d scatter2d = new Scatter2d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, true);
				scatter2d.pack();
				RefineryUtilities.centerFrameOnScreen(scatter2d);
				scatter2d.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
				scatter2d.setSize(1280, 720);
				scatter2d.setVisible(true);
			} else if (3 == problem.getNumberOfObjectives()) {
				new Scatter3d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, true).plot();
			} else {
				final LineBeyend4d lineBeyend4d = new LineBeyend4d("Dimension", "Fitness", problem.getName(), population.writeObjectivesToMatrix(), indicators);
				lineBeyend4d.pack();
				RefineryUtilities.centerFrameOnScreen(lineBeyend4d);
				lineBeyend4d.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
				lineBeyend4d.setSize(1280, 720);
				lineBeyend4d.setVisible(true);
			}
			logger_.info("Total run time is" + endTime + "ms");
			wfghvCalculateRvea wfg = new wfghvCalculateRvea(population, fun);
			double hv = wfg.calculatewfghv();
			logger_.info(problem.getName() + "\nHyperVolume: "
					+ hv + "\nEPSILON    : "
					+ indicators.getEpsilon(population) + "\nGD         : " + indicators.getGD(population) + "\nIGD        : " + indicators.getCEC_IGD(population) + "\nSpread     : " + indicators.getSpread(population));
		}
	} // main
}// AgMOPSO_main