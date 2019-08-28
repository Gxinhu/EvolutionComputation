package jmetal.metaheuristics.MOEAD;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfgHvPlatEMO;
import jmetal.qualityIndicator.hypeHypervolume.HypeHV;
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

public class MOEAD_SBX_main {
	public static void main(String args[]) throws JMException, ClassNotFoundException, IOException {
		int m = 3;
		for (int fun = 6; fun <= 6; fun++) {

			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			int runtimes = 1;
			QualityIndicator indicators = null;
			boolean wfgIs2d = false;
			for (int i = 0; i < runtimes; i++) {
				problem = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgIs2d).getindicator();
				algorithm = new MOEAD_SBX(problem, i);

				if (m == 3) {
					algorithm.setInputParameter("div1", 13);//N=91
					algorithm.setInputParameter("div2", 0);//N=91
					if (problem.getName() == "DTLZ1") {
						algorithm.setInputParameter("maxEvaluations", 500 * 105);
					} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
						algorithm.setInputParameter("maxEvaluations", 100 * 91);
					} else {
						algorithm.setInputParameter("maxEvaluations", 500 * 91);
					}
					algorithm.setInputParameter("swarmSize", 105);
				} else if (m == 5) {
					algorithm.setInputParameter("div1", 6);//N=210
					algorithm.setInputParameter("div2", 0);
					if (problem.getName() == "DTLZ1") {
						algorithm.setInputParameter("maxEvaluations", 200 * 210);
					} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
						algorithm.setInputParameter("maxEvaluations", 80 * 210);
					} else {
						algorithm.setInputParameter("maxEvaluations", 500 * 210);
					}
				} else if (m == 8) {
					algorithm.setInputParameter("div1", 3);//N=156
					algorithm.setInputParameter("div2", 2);
					if (problem.getName() == "DTLZ1") {
						algorithm.setInputParameter("maxEvaluations", 250 * 156);
					} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
						algorithm.setInputParameter("maxEvaluations", 150 * 156);
					} else {
						algorithm.setInputParameter("maxEvaluations", 750 * 156);
					}
				} else if (m == 10) {
					algorithm.setInputParameter("div1", 3);//N=275
					algorithm.setInputParameter("div2", 2);
					if (problem.getName() == "DTLZ1") {
						algorithm.setInputParameter("maxEvaluations", 200 * 275);
					} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
						algorithm.setInputParameter("maxEvaluations", 100 * 275);
					} else {
						algorithm.setInputParameter("maxEvaluations", 750 * 275);
					}
				} else if (m == 15) {
					algorithm.setInputParameter("div1", 2);//N=135
					algorithm.setInputParameter("div2", 1);
					if (problem.getName() == "DTLZ1") {
						algorithm.setInputParameter("maxEvaluations", 500 * 135);
					} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
						algorithm.setInputParameter("maxEvaluations", 250 * 135);
					} else {
						algorithm.setInputParameter("maxEvaluations", 2000 * 135);
					}
				}

				// Mutation and Crossover for Real codification
				HashMap<String, Double> parameters = new HashMap<String, Double>();
				parameters.put("probability", 1.0);
				parameters.put("distributionIndex", 30.0);
				crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

				parameters = new HashMap<String, Double>();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);


				// Add the operators to the algorithm
				algorithm.addOperator("crossover", crossover);
				algorithm.addOperator("mutation", mutation);

				SolutionSet population = algorithm.execute();
				plot(problem, population, indicators);
				wfgHvPlatEMO wfgHvPlatEMO = new wfgHvPlatEMO(population.writeObjectivesToMatrix(), indicators.getTrueParetoFront());
				double time = System.currentTimeMillis();
				double hv = wfgHvPlatEMO.calculatewfghv();
				HypeHV hype = new HypeHV(population.writeObjectivesToMatrix(), indicators.getTrueParetoFront());
				double hvtime = -time + System.currentTimeMillis();
				System.out.println(String.format("The calculate time is %f ", hvtime));
				double hv1 = hype.calculatewfghv();
				double hvtime1 = -time - hvtime + System.currentTimeMillis();
				System.out.println(String.format("The calculate time is %f ", hvtime1));
				assert indicators != null;
				Logger logger = Configuration.getLogger_();
				FileHandler fileHandler = new FileHandler("r2pso.log");
				logger.addHandler(fileHandler);
				logger.info(problem.getName()
						+ "\nHyperVolume: " + hv
						+ "\nHyperVolume1: " + hv1
						+ "\nEPSILON    : " + indicators.getEpsilon(population)
						+ "\nGD         : " + indicators.getGD(population)
						+ "\nIGD        : " + indicators.getCEC_IGD(population)
						+ "\nSpread     : " + indicators.getGeneralizedSpread(population)
						+ "\nSpace        : " + indicators.getSpace(population)
						+ "\nNumberOfPF        : " + population.size()
						+ "\nPD                : " + indicators.getPD(population)
				);
			}

		}//for runtimes
	}//for fun

	public static void plot(Problem problem, SolutionSet population, QualityIndicator indicators) {
		if (2 == problem.getNumberOfObjectives()) {
			final Scatter2d demo = new Scatter2d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, true);
			demo.pack();
			RefineryUtilities.centerFrameOnScreen(demo);
			demo.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
			demo.setSize(1000, 720);
			demo.setVisible(true);
		} else if (3 == problem.getNumberOfObjectives()) {
			new Scatter3d("x", "y", problem.getName(), population.writeObjectivesToMatrix(), indicators, true).plot();
		} else {
			final LineBeyend4d demo = new LineBeyend4d("Dimension", "Fitness", problem.getName(), population.writeObjectivesToMatrix(), indicators);
			demo.pack();
			RefineryUtilities.centerFrameOnScreen(demo);
			demo.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
			demo.setSize(1280, 720);
			demo.setVisible(true);
		}
	}
}//class
