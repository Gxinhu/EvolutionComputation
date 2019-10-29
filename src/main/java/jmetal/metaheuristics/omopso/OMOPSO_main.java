//  OMOPSO_main.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.omopso;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.mutation.Mutation;
import jmetal.operators.mutation.NonUniformMutation;
import jmetal.operators.mutation.UniformMutation;
import jmetal.problems.DTLZ.*;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.WFG.*;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
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

/**
 * Class for configuring and running the OMOPSO algorithm
 */
public class OMOPSO_main {
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
			Configuration.getLogger_().severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGD

	public static void main(String[] args) throws JMException, IOException,
			ClassNotFoundException {
		logger_ = Configuration.getLogger_();
		fileHandler_ = new FileHandler("OMOPSO_main.log");
		logger_.addHandler(fileHandler_);

		for (int fun = 13; fun <= 24; fun++) {
			int runtimes = 30;
			double[] GDarray = new double[runtimes];
			double[] IGDarray = new double[runtimes];
			double[] spreadarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			long Execution_time = 0;
			//long estimatedTime;

			for (int i = 0; i < runtimes; i++) {
				Problem problem = null; // The problem to solve
				Algorithm algorithm; // The algorithm to use
				Mutation uniformMutation;
				Mutation nonUniformMutation;

				QualityIndicator indicators; // Object to get quality indicators

				HashMap parameters; // Operator parameters

				// Logger object and file to store log messages


				indicators = null;
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
					if (fun == 1) {
						problem = new ZDT1("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT1_501.txt");
					}//problem = new WFG1("Real");
					if (fun == 2) {
						problem = new ZDT2("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT2_501.txt");
					}//problem = new WFG1("Real");
					if (fun == 3) {
						problem = new ZDT3("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT3_269.txt");
					}
					if (fun == 4) {
						problem = new ZDT4("Real", 10);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT4_501.txt");
					}//problem = new WFG1("Real");
					if (fun == 5) {
						problem = new ZDT6("Real", 10);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT6_774.txt");
					}//problem = new WFG1("Real");
					if (fun == 6) {
						//problem = new DTLZ1("Real",10,2);
						problem = new DTLZ1("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ1.pf");
					}
					if (fun == 7) {
						//problem = new DTLZ2("Real",10,2);
						problem = new DTLZ2("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ2.pf");
					}//problem = new WFG1("Real");
					if (fun == 8) {
						problem = new DTLZ3("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ3.pf");
					}//problem = new WFG1("Real");
					if (fun == 9) {
						problem = new DTLZ4("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ4.pf");
					}
					if (fun == 10) {
						problem = new DTLZ5("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ5.txt");
					}//problem = new WFG1("Real");
					if (fun == 11) {
						problem = new DTLZ6("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ6.txt");
					}//problem = new WFG1("Real");
					if (fun == 12) {
						problem = new DTLZ7("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ7.pf");
					}
					if (fun == 13) {
						problem = new WFG1("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG1_605.txt");
					}//problem = new WFG1("Real");
					if (fun == 14) {
						problem = new WFG2("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG2_111.txt");
					}//problem = new WFG1("Real");
					if (fun == 15) {
						problem = new WFG3("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG3_right.txt");
					}
					if (fun == 16) {
						problem = new WFG4("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG4_1181.txt");
					}//problem = new WFG1("Real");
					if (fun == 17) {
						problem = new WFG5("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG5_694.txt");
					}
					if (fun == 18) {
						problem = new WFG6("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\wfg6_166.txt");
					}//problem = new WFG1("Real");
					if (fun == 19) {
						problem = new WFG7("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG7_2435.txt");
					}//problem = new WFG1("Real");
					if (fun == 20) {
						problem = new WFG8("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8_201.txt");
					}
					if (fun == 21) {
						problem = new WFG9("Real", 4, 8, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG9_2591.txt");
					}//problem = new WFG1("Real");
					if (fun == 22) {
						problem = new Fonseca("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\Fonseca.pf");
					}//problem = new WFG1("Real");
					if (fun == 23) {
						problem = new Kursawe("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\Kursawe.pf");
					}
					if (fun == 24) {
						problem = new Schaffer("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\Schaffer.pf");
					}//problem = new WFG1("Real");
					if (fun == 25) {
						problem = new UF1("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF1_500.txt");//.txt
					}
					if (fun == 26) {
						problem = new UF2("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF2_500.txt");
					}
					if (fun == 27) {
						problem = new UF3("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF3_500.txt");
					}
					if (fun == 28) {
						problem = new UF4("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF4_500.txt");
					}
					if (fun == 29) {
						problem = new UF5("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF5_21.txt");
					}
					if (fun == 30) {
						problem = new UF6("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF6_668.txt");
					}
					if (fun == 31) {
						problem = new UF7("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF7_500.txt");
					}
					if (fun == 32) {
						problem = new UF8("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF8.DAT");
					}
					if (fun == 33) {
						problem = new UF9("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\UF9.DAT");
					}
				}

				algorithm = new OMOPSO(problem);

				Integer maxIterations = 300;
				Double perturbationIndex = 0.5;
				Double mutationProbability = 1.0 / problem.getNumberOfVariables();

				// Algorithm parameters
				algorithm.setInputParameter("swarmSize", 200);
				algorithm.setInputParameter("archiveSize", 200);
				algorithm.setInputParameter("maxIterations", maxIterations);

				parameters = new HashMap();
				parameters.put("probability", mutationProbability);
				parameters.put("perturbation", perturbationIndex);
				uniformMutation = new UniformMutation(parameters);

				parameters = new HashMap();
				parameters.put("probability", mutationProbability);
				parameters.put("perturbation", perturbationIndex);
				parameters.put("maxIterations", maxIterations);
				nonUniformMutation = new NonUniformMutation(parameters);

				// Add the operators to the algorithm
				algorithm.addOperator("uniformMutation", uniformMutation);
				algorithm.addOperator("nonUniformMutation", nonUniformMutation);

				// Execute the Algorithm
				long initTime = System.currentTimeMillis();
				SolutionSet population = algorithm.execute();
				Execution_time += (System.currentTimeMillis() - initTime);

				// Print the results
		/*logger_.info("Total execution time: " + estimatedTime + "ms");
		logger_.info("Variables values have been writen to file VAR");
		population.printVariablesToFile("VAR");
		logger_.info("Objectives values have been writen to file FUN");
		population.printObjectivesToFile("FUN");

		if (indicators != null) {
			logger_.info("Quality indicators");
			logger_.info("Hypervolume: "
					+ indicators.getHypervolume(population));
			logger_.info("GD         : " + indicators.getGD(population));
			logger_.info("IGD        : " + indicators.getIGD(population));
			logger_.info("Spread     : " + indicators.getSpread(population));
			logger_.info("Epsilon    : " + indicators.getEpsilon(population));
		} // if*/
				//population.printObjectivesToFile("FUN"+fun+"T"+(i+1));
				//population.printObjectivesToFile("Run"+ i + "-FUN-" + problem.getName()
				//		+ "-OMOPSO-FUN");
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
			double sumHypervolume = 0.0;
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
	}// main
} // OMOPSO_main
