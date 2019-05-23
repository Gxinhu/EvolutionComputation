//  NSGAII_main.java
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

package jmetal.metaheuristics.nsgaII;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.DTLZ.*;
import jmetal.problems.ProblemFactory;
import jmetal.problems.WFG.*;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.qualityIndicator.fastHypervolume.wfg.*;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * Class to configure and execute the NSGA-II algorithm.
 * <p>
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba "On the Effect
 * of the Steady-State Selection Scheme in Multi-Objective Genetic Algorithms"
 * 5th International Conference, EMO 2009, pp: 183-197. April 2009)
 */

public class NSGA_DE_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException Usage: three options -
	 *                           jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *                           jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *                           jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *                           paretoFrontFile
	 */
	public static void printGD(String path, double[] GD) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path);//java�ļ�������������ļ���
			OutputStreamWriter osw = new OutputStreamWriter(fos);//OutputStreamWriter���ַ���ͨ���ֽ���������
			BufferedWriter bw = new BufferedWriter(osw);//������
			for (int i = 0; i < GD.length; i++) {
				bw.write(GD[i] + " ");//д��������
				bw.newLine(); //����
			}

			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGD

	public static void printave(String path, double aveIGD, double varianceIGD, double aveHypervolume, double varianceHV) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw = new BufferedWriter(osw);

			// for (int i = 0; i < IGD.length; i++) {

			bw.write(aveIGD + " ");
			bw.newLine();
			bw.write(varianceIGD + " ");
			bw.newLine();
			bw.write(aveHypervolume + " ");
			bw.newLine();
			bw.write(varianceHV + " ");
			bw.newLine();
			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printave

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException, NullPointerException {

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);

		for (int fun = 1; fun <= 31; fun++) {
			int runtimes = 10;
			double[] IGDarray = new double[runtimes];
			long Execution_time = 0;


			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			Operator selection; // Selection operator

			HashMap parameters; // Operator parameters

			QualityIndicator indicators; // Object to get quality indicators
			int m = 3;
			indicators = null;
			boolean wfgis2d = true;
			//choose the problem
			if (args.length == 1) {
				Object[] params = {"Real"};
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = {"Real"};
				indicators = new QualityIndicator(problem, args[1]);
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else { // Default problem
				problem = new cricleselectproblem(problem, indicators, fun, m, wfgis2d).getProblem();
				indicators = new cricleselectproblem(problem, indicators, fun, m, wfgis2d).getindicator();
			}
			for (int i = 0; i < runtimes; i++) {
				algorithm = new NSGA_DE(problem, true, i);
				// Algorithm parameters
				if (problem.getNumberOfObjectives() == 2) {
					if (fun < 6) {
						algorithm.setInputParameter("maxIterations", 250);
					} else if (fun < 22) {
						algorithm.setInputParameter("maxIterations", 500);
					} else {
						algorithm.setInputParameter("maxIterations", 3000);
					}
					algorithm.setInputParameter("swarmSize", 100);
					// Clone operator
				} else if (problem.getNumberOfObjectives() == 3) {
					if (fun < 22) {
						algorithm.setInputParameter("maxIterations", 500);
					} else {
						algorithm.setInputParameter("maxIterations", 3000);
					}
					algorithm.setInputParameter("swarmSize", 105);
					// Clone operator
				}

				algorithm.setInputParameter("T", 10);
				algorithm.setInputParameter("delta", 0.8);
				// Mutation and Crossover for Real codification
				parameters = new HashMap();
				parameters.put("CR", 1.0);
				parameters.put("F", 0.5);
				crossover = CrossoverFactory.getCrossoverOperator(
						"DifferentialEvolutionCrossover", parameters);


				parameters = new HashMap();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation",
						parameters);

				// Selection Operator
				parameters = null;
				selection = SelectionFactory.getSelectionOperator("BinaryTournament3",
						parameters);

				// Add the operators to the algorithm
				algorithm.addOperator("crossover", crossover);
				algorithm.addOperator("mutation", mutation);
				algorithm.addOperator("selection", selection);

				// Add the indicator object to the algorithm
				algorithm.setInputParameter("indicators", indicators);

				// Execute the Algorithm
				long initTime = System.currentTimeMillis();
				SolutionSet population = algorithm.execute();
				Execution_time += (System.currentTimeMillis() - initTime);

				// Result messages
				//population.printVariablesToFile("Variables_NSGAII_3Obj_"+problem.getName()+ "_run"+(i+1)+".txt" );
				IGDarray[i] = indicators.getCEC_IGD(population);
//				wfghvCalculator1 wfg=new wfghvCalculator1(population,fun);
//				double hv=wfg.calculatewfghv();
//				System.out.println(hv);
			}
			double sumIGD = 0;
			for (int i = 0; i < runtimes; i++) {
				sumIGD += IGDarray[i];
			}
			logger_.info("Total execution time: " + Execution_time + "ms");
			// System.out.println("avrHV-fun"+fun+"= "+sumHypervolume/runtimes);
			System.out.println("avrIGD-fun" + fun + "= " + sumIGD / runtimes);
		}
	}

} // NSGAII_main
