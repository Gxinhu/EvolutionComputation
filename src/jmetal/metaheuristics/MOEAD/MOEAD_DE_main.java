package jmetal.metaheuristics.MOEAD;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.DTLZ.*;
import jmetal.problems.ProblemFactory;
import jmetal.problems.WFG.*;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;

public class MOEAD_DE_main {
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

	public static void main(String[] args) throws JMException, ClassNotFoundException {
		int m = 3;
		for (int fun = 1; fun <= 1; fun++) {
			int runtimes = 30;
			if ((fun >= 1 && fun <= 7) || (fun >= 36 && fun <= 44)) {
				m = 3;
			}
			if ((fun >= 8 && fun <= 14) || (fun >= 45 && fun <= 53)) {
				m = 5;
			}
			if ((fun >= 15 && fun <= 21) || (fun >= 54 && fun <= 62)) {
				m = 8;
			}
			if ((fun >= 22 && fun <= 28) || (fun >= 63 && fun <= 71)) {
				m = 10;
			}
			if ((fun >= 29 && fun <= 35) || (fun >= 72 && fun <= 80)) {
				m = 15;
			}
			if (fun == 5 || fun == 12 || fun == 19 || fun == 26 || fun == 33) {
				continue;
			}
			if (fun == 6 || fun == 13 || fun == 20 || fun == 27 || fun == 34) {
				continue;
			}
			if (fun == 7 || fun == 14 || fun == 21 || fun == 28 || fun == 35) {
				continue;
			}
			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator

			if (args.length == 1) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} else if (args.length == 2) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} else { // Default problem
				if (fun == 1 || fun == 8 || fun == 15 || fun == 22 || fun == 29) {
					problem = new DTLZ1("Real", m + 4, m);
				}
				if (fun == 2 || fun == 9 || fun == 16 || fun == 23 || fun == 30) {
					problem = new DTLZ2("Real", m + 9, m);
				}
				if (fun == 3 || fun == 10 || fun == 17 || fun == 24 || fun == 31) {
					problem = new DTLZ3("Real", m + 9, m);
				}
				if (fun == 4 || fun == 11 || fun == 18 || fun == 25 || fun == 32) {
					problem = new DTLZ4("Real", m + 9, m);
				}
				if (fun == 5 || fun == 12 || fun == 19 || fun == 26 || fun == 33) {
					problem = new DTLZ5("Real", m + 9, m);
				}
				if (fun == 6 || fun == 13 || fun == 20 || fun == 27 || fun == 34) {
					problem = new DTLZ6("Real", m + 9, m);
				}
				if (fun == 7 || fun == 14 || fun == 21 || fun == 28 || fun == 35) {
					problem = new DTLZ7("Real", m + 19, m);
				}
				if (fun == 72 || fun == 63 || fun == 54 || fun == 45 || fun == 36) {
					problem = new WFG1("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 73 || fun == 64 || fun == 55 || fun == 46 || fun == 37) {
					problem = new WFG2("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 74 || fun == 65 || fun == 56 || fun == 47 || fun == 38) {
					problem = new WFG3("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 75 || fun == 66 || fun == 57 || fun == 48 || fun == 39) {
					problem = new WFG4("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 76 || fun == 67 || fun == 58 || fun == 49 || fun == 40) {
					problem = new WFG5("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 77 || fun == 68 || fun == 59 || fun == 50 || fun == 41) {
					problem = new WFG6("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 78 || fun == 69 || fun == 60 || fun == 51 || fun == 42) {
					problem = new WFG7("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 79 || fun == 70 || fun == 61 || fun == 52 || fun == 43) {
					problem = new WFG8("Real", 2 * (m - 1), 20, m);
				}
				if (fun == 80 || fun == 71 || fun == 62 || fun == 53 || fun == 44) {
					problem = new WFG9("Real", 2 * (m - 1), 20, m);
				}
			} // else

			algorithm = new MOEAD_DE(problem);

			if (m == 3) {
				algorithm.setInputParameter("div1", 12);//N=91
				algorithm.setInputParameter("div2", 0);//N=91
				if (problem.getName() == "DTLZ1") {
					algorithm.setInputParameter("maxEvaluations", 200 * 91);
				} else if (problem.getName() == "DTLZ2" || problem.getName() == "DTLZ4") {
					algorithm.setInputParameter("maxEvaluations", 100 * 91);
				} else {
					algorithm.setInputParameter("maxEvaluations", 500 * 91);
				}
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
			parameters.put("CR", 0.2);
			parameters.put("F", 0.5);
			crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

			parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			// Add the operators to the algorithm
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);

			for (int i = 0; i < runtimes; i++) {
				SolutionSet population = algorithm.execute();
				population.printObjectivesToFile("results\\MOEAD_DE\\MOEAD_" + problem.getName() + "_" + problem.getNumberOfObjectives() + "_" + "T" + (i + 1));
			}//for runtimes
			System.out.println("problem " + problem.getName() + " with " + problem.getNumberOfObjectives() + " objectives");
		}//for fun
	}//main
}//class
