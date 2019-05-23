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

public class NSGAII_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException,NullPointerException{

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);

		for (int fun = 4; fun <= 4; fun++) {
			int runtimes = 1;
			double[] IGDarray = new double[runtimes];
			long Execution_time = 0;


			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			Operator selection; // Selection operator

			HashMap parameters; // Operator parameters

			QualityIndicator indicators; // Object to get quality indicators
			int m=3;
			indicators=null;
			boolean wfgis2d=true;
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
				problem=new cricleselectproblem(problem,indicators,fun,m,wfgis2d).getProblem();
				indicators=new cricleselectproblem(problem,indicators,fun,m,wfgis2d).getindicator();
			}
			for (int i = 0; i < runtimes; i++) {
				algorithm = new NSGAII(problem,true,i);
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
				parameters.put("probability", 0.9);
				parameters.put("distributionIndex", 20.0);
				crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
						parameters);

				parameters = new HashMap();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation",
						parameters);

				// Selection Operator
				parameters = null;
				selection = SelectionFactory.getSelectionOperator("BinaryTournament",
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
			System.out.println("avrIGD-fun" + fun + "= " + sumIGD / runtimes+problem.getName());
		}
	}
} // NSGAII_main
