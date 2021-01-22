package jmetal.metaheuristics.rvea;
/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgCalRveaExper;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;


public class RveaRuner {
	// Logger object
	public static Logger logger_;
	// FileHandler object
	static FileHandler fileHandler_;

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException {
		// the number of objectives
		int m = 3;
		logger_ = Configuration.getLogger_();
		fileHandler_ = new FileHandler("Rvea.log");
		logger_.addHandler(fileHandler_);
		final int low = 6;
		for (int fun = low; fun <= low; fun++) {
			// The problem to solve
			Problem problem = null;
			// The algorithm to use
			Algorithm algorithm;
			// Crossover operator
			Operator crossover;
			// Mutation operator
			Operator mutation;
			QualityIndicator indicators; // Object to get quality indicators
			indicators = null;
			boolean wfgIs2d = false;
			//choose the problem
			if (args.length == 1) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[
						] params = {"Real"};
				indicators = new QualityIndicator(problem, args[1]);
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else { // Default problem
				problem = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getindicator();
			}
			// init parameter of algorithm
			int i = 0;
			algorithm = new rvea(problem, indicators, i);

			if (fun == 6 | fun == 8) {
				algorithm.setInputParameter("maxIterations", 500);
			} else if (fun <= 12) {
				algorithm.setInputParameter("maxIterations", 500);
			} else if (fun > 12 & fun < 22) {
				algorithm.setInputParameter("maxIterations", 500);
			}
			if (problem.getNumberOfObjectives() == 2) {
				algorithm.setInputParameter("swarmSize", 100);
			} else if (problem.getNumberOfObjectives() == 3) {
				algorithm.setInputParameter("swarmSize", 91);
			} else if (problem.getNumberOfObjectives() == 5) {
				algorithm.setInputParameter("swarmSize", 210);
			} else if (problem.getNumberOfObjectives() == 6) {
				algorithm.setInputParameter("swarmSize", 132);
			} else if (problem.getNumberOfObjectives() == 8) {
				algorithm.setInputParameter("swarmSize", 156);
			} else if (problem.getNumberOfObjectives() == 10) {
				algorithm.setInputParameter("swarmSize", 275);
			}
			HashMap<String, Double> parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0);
			parameters.put("distributionIndex", 30.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

			parameters = new HashMap<>();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			// Add the operators to the algorithm
			SolutionSet population = null;
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			long initTime = System.currentTimeMillis();
			population = algorithm.execute();
			long endTime = System.currentTimeMillis() - initTime;
			//画图
			logger_.info("Total run time is" + endTime + "ms");
			wfgCalRveaExper wfg = new wfgCalRveaExper(population.writeObjectivesToMatrix(), problem.getName());
			double hv = wfg.calculatewfghv();
			logger_.info(problem.getName() + "\nHyperVolume: "
					+ hv + "\nEPSILON    : "
					+ indicators.getEpsilon(population)
					+ "\nGD         : " + indicators.getGD(population)
					+ "\nIGD        : " + indicators.getCEC_IGD(population)
					+ "\nSpace      : " + indicators.getSpace(population)
					+ "\nSpread     : " + indicators.getGeneralizedSpread(population)
					+ "\nPD         : " + indicators.getPD(population));


		}

	} // main
}// AgMOPSO_main