/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */
package jmetal.metaheuristics.AgMOPSO;

import jmetal.core.*;
import jmetal.operators.clone.CloneFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.WFG.*;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF10;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.R2;
import jmetal.qualityIndicator.fastHypervolume.wfg.Front;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.HashMap;


public class AgMOPSO_main {
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
		long startime=System.currentTimeMillis();
		int m = 3;
		for (int fun = 1; fun <=1; fun++) {
			int runtimes = 1;
			double[] IGDarray = new double[runtimes];
			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator clone = null; // Crossover operator
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			QualityIndicator indicators; // Object to get quality indicators
			indicators = null;
			for (int i = 0; i < runtimes; i++) {

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

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/ZDT/ZDT1.pf");
					}//problem = new WFG1("Real");
					if (fun == 2) {
						problem = new ZDT2("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/ZDT/ZDT2.pf");
					}//problem = new WFG1("Real");
					if (fun == 3) {
						problem = new ZDT3("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/ZDT/ZDT3.pf");
					}
					if (fun == 4) {
						problem = new ZDT4("Real", 10);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/ZDT/ZDT4.pf");
					}//problem = new WFG1("Real");
					if (fun == 5) {
						problem = new ZDT6("Real", 10);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/ZDT/ZDT6.pf");
					}//problem = new WFG1("Real");
					if (fun == 6) {
						problem = new DTLZ1("Real", m + 4, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ1.pf");
					}
					if (fun == 7) {
						problem = new DTLZ2("Real", m + 9, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ2.pf");
					}//problem = new WFG1("Real");
					if (fun == 8) {
						problem = new DTLZ3("Real", m + 9, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ3.pf");
					}//problem = new WFG1("Real");
					if (fun == 9) {
						problem = new DTLZ4("Real", m + 9, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ4.pf");
					}
					if (fun == 10) {
						problem = new DTLZ5("Real", m + 9, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ5.pf");
					}//problem = new WFG1("Real");
					if (fun == 11) {
						problem = new DTLZ6("Real", m + 9, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ6.pf");
					}//problem = new WFG1("Real");
					if (fun == 12) {
						problem = new DTLZ7("Real", m + 19, m);

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/DTLZ/" + m + "d/DTLZ7.pf");
					}
//					if (fun == 13) {
//						problem = new WFG1("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG1.pf");
//					}//problem = new WFG1("Real");
//					if (fun == 14) {
//						problem = new WFG2("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG2.pf");
//					}//problem = new WFG1("Real");
//					if (fun == 15) {
//						problem = new WFG3("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG3.pf");
//					}
//					if (fun == 16) {
//						problem = new WFG4("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG4.pf");
//					}//problem = new WFG1("Real");
//					if (fun == 17) {
//						problem = new WFG5("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG5.pf");
//					}
//					if (fun == 18) {
//						problem = new WFG6("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG6.pf");
//					}//problem = new WFG1("Real");
//					if (fun == 19) {
//						problem = new WFG7("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG7.pf");
//					}//problem = new WFG1("Real");
//					if (fun == 20) {
//						problem = new WFG8("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG8.pf");
//					}
//					if (fun == 21) {
//						problem = new WFG9("Real", m - 1, 10, m);
//
//						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG9.pf");
//					}
                    if (fun == 13) {
                        problem = new WFG1("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG1.pf");
                    }//problem = new WFG1("Real");
                    if (fun == 14) {
                        problem = new WFG2("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG2.pf");
                    }//problem = new WFG1("Real");
                    if (fun == 15) {
                        problem = new WFG3("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG3.pf");
                    }
                    if (fun == 16) {
                        problem = new WFG4("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG4.pf");
                    }//problem = new WFG1("Real");
                    if (fun == 17) {
                        problem = new WFG5("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG5.pf");
                    }
                    if (fun == 18) {
                        problem = new WFG6("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG6.pf");
                    }//problem = new WFG1("Real");
                    if (fun == 19) {
                        problem = new WFG7("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG7.pf");
                    }//problem = new WFG1("Real");
                    if (fun == 20) {
                        problem = new WFG8("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG8.pf");
                    }
                    if (fun == 21) {
                        problem = new WFG9("Real", 4, 20, 2);

                        indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/WFG/" + m + "d/WFG9.pf");
                    }
					if (fun == 22) {
						problem = new UF1("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF1.pf");//.txt
					}
					if (fun == 23) {
						problem = new UF2("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF2.pf");
					}
					if (fun == 24) {
						problem = new UF3("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF3.pf");
					}
					if (fun == 25) {
						problem = new UF4("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF4.pf");
					}
					if (fun == 26) {
						problem = new UF5("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF5.pf");
					}
					if (fun == 27) {
						problem = new UF6("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF6.pf");
					}
					if (fun == 28) {
						problem = new UF7("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF7.pf");
					}
					if (fun == 29) {
						problem = new UF8("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF8.pf");
					}
					if (fun == 30) {
						problem = new UF9("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF9.pf");
					}
					if (fun == 31) {
						problem = new UF10("Real");

						indicators = new QualityIndicator(problem, "/home/hu/PSO/jmetal Agmopso/PF/UF/UF10.pf");
					}
				}

				SolutionSet population = null;
//                algorithm = new AgmopsoAdaptiveWeighter(problem, indicators, i);
//                algorithm=new AgR2ADW(problem,indicators,i);
                algorithm=new AgMOPSO(problem,indicators,i);
//				algorithm = new AgMOPSOwithR2newVersion(problem);
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
				parameters.put("probability", 1.0);
				parameters.put("distributionIndex", 20.0);
				crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

				parameters = new HashMap<String, Double>();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

				// Add the operators to the algorithm
				algorithm.addOperator("clone", clone);
				algorithm.addOperator("crossover", crossover);
				algorithm.addOperator("mutation", mutation);

				population = algorithm.execute();
				IGDarray[i] = indicators.getCEC_IGD(population);
			}
			Arrays.sort(IGDarray);
			TestStatistics sta = null;
			sta = new TestStatistics(IGDarray);
			long endtime=System.currentTimeMillis();
			System.out.println(endtime-startime+"ms");
			System.out.println(sta.getAverage() + "\t" + sta.getStandardDiviation() + "\t" + problem.getNumberOfObjectives() + problem.getName());
			//System.out.println(problem.getName());
		} //runtimes
	} // main
}// AgMOPSO_main