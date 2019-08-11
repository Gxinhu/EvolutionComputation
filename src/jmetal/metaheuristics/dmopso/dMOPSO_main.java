/**
 * dMOPSO_main.java
 *
 * @version 1.0
 */

package jmetal.metaheuristics.dmopso;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.mutation.Mutation;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

public class dMOPSO_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException Usage: three options -
	 *                           jmetal.metaheuristics.mocell.MOCell_main -
	 *                           jmetal.metaheuristics.mocell.MOCell_main problemName -
	 *                           jmetal.metaheuristics.mocell.MOCell_main problemName
	 *                           ParetoFrontFile
	 */
	public static void printGD(String path, double[] GD) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);
			for (int i = 0; i < GD.length; i++) {
				bw.write(GD[i] + " ");
				bw.newLine();
			}

			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGD

	public static void main(String[] args) throws JMException, IOException,
			ClassNotFoundException {


		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("dMOPSO_main.log");
		logger_.addHandler(fileHandler_);
		for (int fun = 6; fun <= 6; fun++) {
			int runtimes = 1;
			double[] GDarray = new double[runtimes];
			double[] IGDarray = new double[runtimes];
			double[] spreadarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			long Execution_time = 0;

			for (int i = 0; i < runtimes; i++) {

				Problem problem; // The problem to solve
				Algorithm algorithm; // The algorithm to use
				Mutation mutation; // "Turbulence" operator
				HashMap parameters; // Operator parameters
				QualityIndicator indicators; // Object to get quality indicators
				indicators = null;
				problem = null;
				int m = 3;
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

				} // else

				algorithm = new dMOPSO(problem);

				// Algorithm parameters
				algorithm.setInputParameter("swarmSize", 105);
				algorithm.setInputParameter("maxAge", 2);
				algorithm.setInputParameter("maxIterations", 500);
				algorithm.setInputParameter("functionType", "_PBI");

				parameters = new HashMap();

				// Execute the Algorithm
				long initTime = System.currentTimeMillis();
				SolutionSet population = algorithm.execute();
				Execution_time += (System.currentTimeMillis() - initTime);

				// Result messages
		/*logger_.info("Total execution time: " + estimatedTime + "ms");
		logger_.info("Objectives values have been writen to file FUN");
		population.printObjectivesToFile("FUN");
		logger_.info("Variables values have been writen to file VAR");
		population.printVariablesToFile("VAR");

		if (indicators != null) {
			logger_.info("Quality indicators");
			logger_.info("Hypervolume: "
					+ indicators.getHypervolume(population));
			logger_.info("GD         : " + indicators.getGD(population));
			logger_.info("IGD        : " + indicators.getIGD(population));
			logger_.info("Spread     : " + indicators.getSpread(population));
			logger_.info("Epsilon    : " + indicators.getEpsilon(population));
		} // if*/
				//population.printObjectivesToFile("Run"+ i + "-FUN-" + problem.getName()
				//		+ "-DMOPSO-FUN");
				GDarray[i] = indicators.getCEC_IGD(population);
				IGDarray[i] = indicators.getIGD(population);
				spreadarray[i] = indicators.getEpsilon(population);
				Hypervolume[i] = indicators.getHypervolume(population);
			}
			printGD("FUN" + fun + "GD", GDarray);
			//printGD("FUN"+fun+"IGD",IGDarray);
			//printGD("FUN"+fun+"IGD",IGDarray);
			//printDiversity("FUN"+fun+"Diversity",spreadarray);
			//printHypervolume("FUN"+fun+"Hypervolume",Hypervolume);
			double sumGD = 0;
			double sumIGD = 0;
			double sumSP = 0;
			double sumHypervolume = 0;
			for (int i = 0; i < runtimes; i++) {
				sumGD += GDarray[i];
				sumIGD += IGDarray[i];
				sumSP += spreadarray[i];
				sumHypervolume += Hypervolume[i];
			}
			logger_.info("Total execution time: " + Execution_time + "ms");
			System.out.println(sumGD / runtimes);
			System.out.println(sumSP / runtimes);
			System.out.println(sumIGD / runtimes);
			System.out.println(sumHypervolume / runtimes);
		}
	} // main
} // dMOPSO_main