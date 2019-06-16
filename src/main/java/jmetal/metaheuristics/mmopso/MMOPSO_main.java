package jmetal.metaheuristics.mmopso;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
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

public class MMOPSO_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

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

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		QualityIndicator indicators; // Object to get quality indicators
		Operator mutation; // Mutation operator
		Operator crossover; // Crossover operator
		HashMap parameters; // Operator parameters
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("DDMOPSO_main.log");
		logger_.addHandler(fileHandler_);
		for (int fun = 13; fun <= 21; fun++) {

			problem = null;
			indicators = null;

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
			//problem = new WFG8 ("Real");

			String path = System.getProperty("user.dir") + "/RESULTS-DDMOPSO-"
					+ problem.getName();
			String curDir = System.getProperty("user.dir");
			String str = curDir + "/RESULTS-DDMOPSO-" + problem.getName();


			//createFolder(str);
			//CrowdingArchive leaders_ = new CrowdingArchive(20000,
			//		problem.getNumberOfObjectives());

			SolutionSet population = null;

			//CrowdingArchive finall = new CrowdingArchive(20000,
			//		problem.getNumberOfObjectives());
			int runtimes = 30;
			double[] IGDarray = new double[runtimes];
			double[] spreadarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			long estimatedTime = 0;
			for (int run = 0; run < runtimes; run++) {

				//leaders_ = new CrowdingArchive(10000,
				//		problem.getNumberOfObjectives());

				algorithm = new MMOPSO2(problem);
				//algorithm = new DDMOPSO(problem, run);

				//indicators = new QualityIndicator(problem,"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8.2D.pf" ) ;
				algorithm.setInputParameter("functionType", "_NBI");//_WSUM _NBI
				algorithm.setInputParameter("swarmSize", 200);
				algorithm.setInputParameter("maxIterations", 300);
				algorithm.setInputParameter("selectionArchiveSize", 200
				);
				parameters = new HashMap();
				parameters.put("probability", 1.0);
				parameters.put("distributionIndex", 20.0);
				crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
						parameters);

				parameters = new HashMap();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation",
						parameters);

				algorithm.addOperator("mutation", mutation);

				algorithm.addOperator("crossover", crossover);

				long initTime = System.currentTimeMillis();
				population = algorithm.execute();
				estimatedTime += System.currentTimeMillis() - initTime;

				for (int i = 0; i < population.size(); i++) {
					//leaders_.add(population.get(i));
					//finall.add(population.get(i));
				}
				//System.out.println("RUN-" + run + " is DONE");
				IGDarray[run] = indicators.getCEC_IGD(population);
				spreadarray[run] = indicators.getEpsilon(population);
				Hypervolume[run] = indicators.getHypervolume(population);
				//population.printObjectivesToFile("Run"+ run + "-FUN-" + problem.getName()
				//		+ "-MSMPSO-FUN");
				//finall.printObjectivesToFile("FINAL-FUN-" + problem.getName()
				//		+ "-DDMOPSO-FUN");
			}
			//long estimatedTime = System.currentTimeMillis() - initTime;
			// Result messages
			double sumIGD = 0;
			double sumSP = 0;
			double sumHypervolume = 0.0;
			for (int i = 0; i < runtimes; i++) {
				sumIGD += IGDarray[i];
				sumSP += spreadarray[i];
				sumHypervolume += Hypervolume[i];
			}
			System.out.println(sumIGD / runtimes);
			System.out.println(sumSP / runtimes);
			System.out.println(sumHypervolume / runtimes);
			printGD("FUN" + fun + "IGD", IGDarray);
			logger_.info("Total execution time: " + estimatedTime + "ms");
			//logger_.info("Objectives values have been writen to file FUN");

		}

	} // main

} // MOPSOD_main