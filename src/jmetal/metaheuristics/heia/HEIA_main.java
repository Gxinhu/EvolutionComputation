package jmetal.metaheuristics.heia;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
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
 * Class to configure and execute the NCIA algorithm.
 */

public class HEIA_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException Usage: three options - jmetal.metaheuristics.nica.NICA_main -
	 *                           jmetal.metaheuristics.nica.NICA_main problemName -
	 *                           jmetal.metaheuristics.nica.NICA_main problemName
	 *                           paretoFrontFile
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

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NICA_main.log");
		logger_.addHandler(fileHandler_);


		int fun;
		for (fun = 13; fun <= 14; fun++) {
			int runtimes = 30;
			double[] GDarray = new double[runtimes];
			double[] IGDarray = new double[runtimes];
			double[] spreadarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			long Execution_time = 0;


			for (int i = 0; i < runtimes; i++) {
				Problem problem = null; // The problem to solve
				Algorithm algorithm; // The algorithm to use
				Operator clone = null; // Crossover operator
				Operator crossover, crossover2; // Crossover operator
				Operator mutation; // Mutation operator
				//Operator selection ; // Selection operator

				HashMap parameters; // Operator parameters

				QualityIndicator indicators; // Object to get quality indicators

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
					// problem = new Kursawe("Real", 3);
					// problem = new ZDT1("ArrayReal", 100);
					// problem = new Kursawe("BinaryReal", 3);
					// problem = new Water("Real");
					//problem = new UF1("Real", 30);
					//indicators = new QualityIndicator(problem,
					//"E:\\new_multiobjective\\jMetal\\Pareto_front\\UF1.DAT");
					// problem = new ConstrEx("Real");
					// problem = new DTLZ1("Real");
					// problem = new OKA2("Real") ;
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
						problem = new DTLZ1("Real", 10, 3);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ1.pf");
					}
					if (fun == 7) {
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
						problem = new DTLZ5("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ5.pf");
					}//problem = new WFG1("Real");
					if (fun == 11) {
						problem = new DTLZ6("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ6.pf");
					}//problem = new WFG1("Real");
					if (fun == 12) {
						problem = new DTLZ7("Real");

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ7.pf");
					}
					if (fun == 13) {
						problem = new WFG1("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "/home/hu/Desktop/EvolutionComputation/PF/WFG/2d/WFG1.pf");
					}//problem = new WFG1("Real");
					if (fun == 14) {
						problem = new WFG2("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "/home/hu/Desktop/EvolutionComputation/PF/WFG/2d/WFG2.pf");
					}//problem = new WFG1("Real");
					if (fun == 15) {
						problem = new WFG3("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG3_right.txt");
					}
					if (fun == 16) {
						problem = new WFG4("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG4_1181.txt");
					}//problem = new WFG1("Real");
					if (fun == 17) {
						problem = new WFG5("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG5_694.txt");
					}
					if (fun == 18) {
						problem = new WFG6("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\wfg6_166.txt");
					}//problem = new WFG1("Real");
					if (fun == 19) {
						problem = new WFG7("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG7_2435.txt");
					}//problem = new WFG1("Real");
					if (fun == 20) {
						problem = new WFG8("Real", 8, 2, 2);

						indicators = new QualityIndicator(problem, "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8.2D.pf");
					}
					if (fun == 21) {
						problem = new WFG9("Real", 8, 2, 2);

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
				} // else

				algorithm = new HEIA(problem);

				// Algorithm parameters

				if (fun >= 1 && fun <= 5) {
					algorithm.setInputParameter("populationSize", 100);
					algorithm.setInputParameter("maxEvaluations", 25000);
					// Clone operator
					parameters = new HashMap();
					parameters.put("clonesize", 100);
					clone = CloneFactory.getClone("proportionalclone", parameters);
					//clone = CloneFactory.getClone("proportional2", parameters);
				} else if (fun >= 6 && fun <= 12) {
					algorithm.setInputParameter("populationSize", 500);
					algorithm.setInputParameter("maxEvaluations", 100000);
					// Clone operator
					parameters = new HashMap();
					parameters.put("clonesize", 500);
					clone = CloneFactory.getClone("proportionalclone", parameters);
					//clone = CloneFactory.getClone("proportional2", parameters);
				} else if (fun >= 13 && fun <= 21) {
					algorithm.setInputParameter("populationSize", 200);
					algorithm.setInputParameter("maxEvaluations", 100000);
					// Clone operator
					parameters = new HashMap();
					parameters.put("clonesize", 100);
					clone = CloneFactory.getClone("proportionalclone", parameters);
					//clone = CloneFactory.getClone("proportional2", parameters);
				} else if (fun >= 25 && fun <= 31) {
					algorithm.setInputParameter("populationSize", 300);
					algorithm.setInputParameter("maxEvaluations", 300000);
					// Clone operator
					parameters = new HashMap();
					parameters.put("clonesize", 300);
					clone = CloneFactory.getClone("proportionalclone", parameters);
					//clone = CloneFactory.getClone("proportional2", parameters);
				}

	/*
		// Clone operator
		parameters = new HashMap();
		parameters.put("clonesize", 200);
		clone = CloneFactory.getClone("proportionalclone", parameters);
		//clone = CloneFactory.getClone("proportional2", parameters);
	*/
				// Mutation and Crossover for Real codification

				parameters = new HashMap();
				parameters.put("CR", 1.0);
				parameters.put("F", 0.5);
				crossover = CrossoverFactory.getCrossoverOperator(
						"DifferentialEvolutionCrossover", parameters);


				parameters = new HashMap();
				parameters.put("probability", 1.0);
				parameters.put("distributionIndex", 20.0);
				crossover2 = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
		 
		/*parameters = new HashMap();
		parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		mutation = MutationFactory.getMutationOperator("AdaptiveMutation",
				parameters);*/
				parameters = new HashMap();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation",
						parameters);

				// Selection Operator

				// parameters = null ;
				//selection =SelectionFactory.getSelectionOperator("BinaryTournament2",
				// parameters) ;


				// Add the operators to the algorithm
				algorithm.addOperator("clone", clone);
				algorithm.addOperator("crossover", crossover);
				algorithm.addOperator("crossover2", crossover2);
				algorithm.addOperator("mutation", mutation);
				//algorithm.addOperator("selection",selection);

				// Add the indicator object to the algorithm
				algorithm.setInputParameter("indicators", indicators);

				// Execute the Algorithm
				long initTime = System.currentTimeMillis();
				SolutionSet population = algorithm.execute();
				//long estimatedTime = System.currentTimeMillis() - initTime;

				Execution_time += (System.currentTimeMillis() - initTime);

				//logger_.info("Objectives values have been writen to file FUN"+(i+1));
				//population.printObjectivesToFile("FUN"+(i+1));
				// Result messages
//				logger_.info("Total execution time: " + Execution_time+ "ms");
		/*logger_.info("Variables values have been writen to file VAR");
		// population.printVariablesToFile("VAR");
		logger_.info("Objectives values have been writen to file FUN");
		population.printObjectivesToFile("FUN");

		if (indicators != null) {
			logger_.info("Quality indicators");
			logger_.info("Hypervolume: "
					+ indicators.getHypervolume(population));
			logger_.info("GD         : " + indicators.getGD(population));
			logger_.info("CECIGD        : " + indicators.getCECIGD(population));
			logger_.info("Spread     : " + indicators.getSpread(population));
			// logger_.info("Epsilon    : " + indicators.getEpsilon(population))
			// ;

			int evaluations = ((Integer) algorithm
					.getOutputParameter("evaluations")).intValue();
			logger_.info("Speed      : " + evaluations + " evaluations");
		} // if*/
				//population.printObjectivesToFile("FUN"+fun+"T"+(i+1) + "-HEIA");
				//GDarray[i]=indicators.getGD(population);
				IGDarray[i] = indicators.getCEC_IGD(population);
				//spreadarray[i]=indicators.getSpread(population);
				//Hypervolume[i]=indicators.getHypervolume(population);
			}
			//printGD("FUN"+fun+"GD",GDarray);
			printGD("FUN" + fun + "IGD" + "-HEIA", IGDarray);
			//printDiversity("FUN"+fun+"Diversity",spreadarray);
			//printHypervolume("FUN"+fun+"Hypervolume",Hypervolume);
			double sumGD = 0;
			double sumIGD = 0;
			double sumSP = 0;
			for (int i = 0; i < runtimes; i++) {
				sumGD += GDarray[i];
				sumIGD += IGDarray[i];
				sumSP += spreadarray[i];
			}
			//logger_.info("Total execution time: " + Execution_time + "ms");
			//System.out.println(sumGD/runtimes);
			//System.out.println(sumSP/runtimes);
			System.out.println(sumIGD / runtimes);
		}
	} // main
} // NICA_main