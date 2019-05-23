/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */
package jmetal.metaheuristics.AgMOPSO;

import jmetal.core.*;
import jmetal.metaheuristics.cricleselectproblem;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.qualityIndicator.fastHypervolume.wfg.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;


public class Agmopsorunner {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException,NullPointerException{
		// the numbes of objectives
		int m = 3;
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("MOEAD.log");
		logger_.addHandler(fileHandler_);
		for (int fun = 1; fun <= 1; fun++) {
			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator clone = null; // Crossover operator
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			QualityIndicator indicators; // Object to get quality indicators
			indicators = null;
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
			// init parameter of algorithm
			int i = 0;
//                algorithm = new AgmopsoAdaptiveWeighter(problem, indicators, i);
//                algorithm=new AgR2ADW(problem,indicators,i);
			algorithm = new AgmopsoDE(problem, indicators, i,false);
//				algorithm = new AgMOPSOwithR2newVersion(problem);
//			algorithm = new AgMOPSO(problem, indicators, i,false);


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
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 100);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 3) {
				if (fun < 22) {
					algorithm.setInputParameter("maxIterations", 500);
				} else {
					algorithm.setInputParameter("maxIterations", 3000);
				}
				algorithm.setInputParameter("swarmSize", 105);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 105);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 5) {
				algorithm.setInputParameter("maxIterations", 500);
				algorithm.setInputParameter("swarmSize", 210);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 210);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 8) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 156);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 156);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			} else if (problem.getNumberOfObjectives() == 10) {
				algorithm.setInputParameter("maxIterations", 1000);
				algorithm.setInputParameter("swarmSize", 275);
				// Clone operator
				HashMap<String, Integer> parameters = new HashMap<String, Integer>();
				parameters.put("clonesize", 275);
				clone = CloneFactory.getClone("proportionalclone", parameters);
			}
			HashMap<String, Double> parameters = new HashMap<String, Double>();
//			parameters.put("probability", 1.0);
//			parameters.put("distributionIndex", 20.0);
//			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
			parameters.put("CR", 1.0);
			parameters.put("F", 0.5);
			crossover = CrossoverFactory.getCrossoverOperator(
					"DifferentialEvolutionCrossover", parameters);

			parameters = new HashMap<String, Double>();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

			// Add the operators to the algorithm
			SolutionSet population = null;
			algorithm.addOperator("clone", clone);
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			long initime=System.currentTimeMillis();
			population = algorithm.execute();
			long endtime=System.currentTimeMillis()-initime;
			logger_.info("Total run time is"+endtime+"ms");
			wfghvCalculator1 wfg=new wfghvCalculator1(population);
			double hv=wfg.calculatewfghv();
			if (indicators != null) {
				logger_.info(problem.getName()+"\nHypervolume: "
						+ hv+"\nEPSILON    : "
						+ indicators.getEpsilon(population)+"\nGD         : " + indicators.getGD(population)+"\nIGD        : " + indicators.getCEC_IGD(population)+"\nSpread     : " + indicators.getSpread(population));
			}
		}
//runtimes
	} // main
}// AgMOPSO_main