/**
 * VaPSORunner.java
 *
 * @author Xin.Hu
 */
package jmetal.dynamicAlogrithms.CMPSODMO;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.dynamicAlogrithms.selectproblem;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgHvPlatEMO;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.plot.pythonplot;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;


public class CMPSODMORunner {

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException, InterruptedException {
		// the numbers of objectives
		int m = 3;
		final int low = 1;
		final int high = 1;
		final int yitat = 10;
		final int taut = 10;
		final int t0 = 10;
		Logger logger = Configuration.getLogger_();
		FileHandler fileHandler = new FileHandler("Vepso.log");
		logger.addHandler(fileHandler);
		for (int fun = low; fun <= high; fun++) {
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
				problem = new selectproblem(problem, fun, m, yitat, taut, t0).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getindicator();

			}
			// init parameter of algorithm
			algorithm = new CMPSODMO(problem);
			coffientSetting(algorithm, problem, fun);
			SolutionSet population;
			long initTime = System.currentTimeMillis();
			population = algorithm.execute();
			long endTime = System.currentTimeMillis() - initTime;
			logger.info("Total run time is" + endTime + "ms");

			wfgHvPlatEMO wfgHvPlatEMO = new wfgHvPlatEMO(population.writeObjectivesToMatrix(), problem.getName());
			double hv = wfgHvPlatEMO.calculatewfghv();
			System.out.println(population.size());
			pythonplot plot = new pythonplot(population.writeObjectivesToMatrix(), problem.getName());
			plot.exectue();
//			plot = new pythonplot(problem.getPF(), problem.getName());
//			plot.exectue();
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
				algorithm.setInputParameter("maxIterations", 100);
			} else if (fun < 22) {
				algorithm.setInputParameter("maxIterations", 150);
			} else {
				algorithm.setInputParameter("maxIterations", 3000);
			}
			algorithm.setInputParameter("swarmSizes", 20);
		} else if (problem.getNumberOfObjectives() == 3) {
			if (fun < 22) {
				algorithm.setInputParameter("maxIterations", 99);
			} else {
				algorithm.setInputParameter("maxIterations", 3000);
			}
			algorithm.setInputParameter("swarmSizes", 20);
		}
		algorithm.setInputParameter("archiveSize", 100);

	}

}