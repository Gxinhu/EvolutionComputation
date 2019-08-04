/**
 * NSGAIIStudy.java
 * This is the latest experimental study of NSGAII
 */

package jmetal.experiments.studies;

import jmetal.core.Algorithm;
import jmetal.experiments.Experiment;
import jmetal.experiments.Settings;
import jmetal.experiments.settings.r2psoSetting;
import jmetal.experiments.util.Friedman;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;


public class myStudy extends Experiment {

	/**
	 * Configures the algorithms in each independent run
	 *
	 * @param problemName  The problem to solve
	 * @param problemIndex The index of the problem
	 */
	@Override
	public void algorithmSettings(String problemName,
	                              int problemIndex,
	                              Algorithm[] algorithm) throws ClassNotFoundException {
		try {
			int numberOfAlgorithms = algorithmNameList_.length;

			HashMap[] parameters = new HashMap[numberOfAlgorithms];

			Object[] problemParams;

			for (int i = 0; i < numberOfAlgorithms; i++) {
				parameters[i] = new HashMap();
			} // for

			if (!paretoFrontFile_[problemIndex].equals("")) {
				for (int i = 0; i < numberOfAlgorithms; i++) {
					parameters[i].put("paretoFrontFile_", paretoFrontFile_[problemIndex]);
				}
			} // if
			/**
			 * Real code
			 */
			if (problemName.contains("WFG")) {
				// WFG: solutionType, k, l, M
				int k = 0;
				k = 2 * (noOfObjectives_ - 1);
				problemParams = new Object[]{"Real", k, 20, noOfObjectives_};


			} else if (problemName.contains("DTLZ")) {
				int n;
				if (problemName.equalsIgnoreCase("DTLZ1") || problemName.equalsIgnoreCase("ScaledDTLZ1")) {
					// r = 5 for DTLZ1 serious
					n = noOfObjectives_ + 5 - 1;
				} else if (problemName.equalsIgnoreCase("DTLZ7")) {
					// r = 20 for DTLZ7
					n = noOfObjectives_ + 20 - 1;
				} else {
					// r = 10 for DTLZ2-6
					n = noOfObjectives_ + 10 - 1;
				}
				// DTLZ: solutionType,  n, m
				problemParams = new Object[]{"Real", n, noOfObjectives_};

			} else if (problemName.contains("MaF")) {
				int n;
				if (problemName.equalsIgnoreCase("MaF2")
						|| problemName.equalsIgnoreCase("MaF3") || problemName.equalsIgnoreCase("MaF4")
						|| problemName.equalsIgnoreCase("MaF5") || problemName.equalsIgnoreCase("MaF6")) {

					// D = M + K - 1, K = 10
					n = noOfObjectives_ + 10 - 1;

				} else if (problemName.equalsIgnoreCase("MaF1")) {
					n = noOfObjectives_ + 5 - 1;
				} else if (problemName.equalsIgnoreCase("MaF7")) {
					// K = 20 for MaF7
					n = noOfObjectives_ + 20 - 1;

				} else if (problemName.equalsIgnoreCase("MaF8") || problemName.equalsIgnoreCase("MaF9")) {
					// MaF8 and MaF9  D = 2;
					n = 2;

				} else if (problemName.equalsIgnoreCase("MaF10") || problemName.equalsIgnoreCase("MaF11")
						|| problemName.equalsIgnoreCase("MaF12")) {
					// MaF10--12  L=10, D = M - 1 + L
					n = noOfObjectives_ - 1 + 10;

				} else if (problemName.equalsIgnoreCase("MaF13")) {
					// MaF13 D = 5
					n = 5;

				} else {
					// MaF14 and 15 D = 20 * M
					n = 20 * noOfObjectives_;
				}

				// solutionType,  n, m
				problemParams = new Object[]{"Real", n, noOfObjectives_};

			} else {

				problemParams = new Object[]{"Real"};

			}
//			algorithm[0] = new RVEASettting(problemName, problemParams).configure(parameters[0]);
			algorithm[0] = new r2psoSetting(problemName, problemParams).configure(parameters[0]);
//			algorithm[1] = new NSGAIIISetting(problemName, problemParams).configure(parameters[0]);
			//			algorithm[1] = new r2psoSetting(problemName, problemParams).configure(parameters[0]);

		} catch (IllegalArgumentException | IllegalAccessException | JMException ex) {
			Logger.getLogger(RVEAStudy.class.getName()).log(Level.SEVERE, null, ex);
		}
	} // algorithmSettings

	/**
	 * Main method
	 *
	 * @param args
	 * @throws JMException
	 * @throws IOException
	 */
	public static void main(String[] args) throws JMException, IOException {
		int[] index = {3};
		for (int j = 0; j < index.length; j++) {
			myStudy exp = new myStudy();

			exp.experimentName_ = "VagPSOStudy";

			exp.algorithmNameList_ = new String[]{
					"VagPSOorigin",
					"VagPSO_addBigAngle",
//					"VagPSOaddtheBestnorm",
					"VagPSOclonetheEDge",
					"VagPSO"
			};
			exp.problemList_ = new String[]{
					"DTLZ1",
					"DTLZ2",
					"DTLZ3",
					"DTLZ4",
					"DTLZ5",
					"DTLZ6",
					"DTLZ7",
					"WFG1",
					"WFG2",
					"WFG3",
					"WFG4",
					"WFG5",
					"WFG6",
					"WFG7",
					"WFG8",
					"WFG9",
//				"MaF1",
//				"MaF2",
//				"MaF3",
//				"MaF4",
//				"MaF5",
//				"MaF6",
//				"MaF7",
					/**
					 * M = 2
					*/
//			"UF1","UF2",//"UF3","UF4","UF5","UF6","UF7",
//			//"ZDT1","ZDT2","ZDT3",
//			"ZDT4","ZDT6",
			};

			exp.noOfObjectives_ = index[j];
			exp.paretoFrontFile_ = new String[]{
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ1.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ2.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ3.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ4.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ5.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ6.pf",
					"DTLZ/" + exp.noOfObjectives_ + "d/DTLZ7.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG1.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG2.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG3.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG4.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG5.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG6.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG7.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG8.pf",
					"WFG/" + exp.noOfObjectives_ + "d/WFG9.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF1.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF2.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF3.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF4.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF5.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF6.pf",
//				"Maf/" + exp.noOfObjectives_ + "d/MaF7.pf",
			};

			exp.indicatorList_ = new String[]{

//				"HV",
//				"HV2",
//				"GSPREAD",
//    		"DCI",
//    		"EPSILON",
					"IGD",
//				"Space",
//				"PD",
//				"GD",
//    		"RUNTIME",

			};

			int numberOfAlgorithms = exp.algorithmNameList_.length;

			exp.experimentBaseDirectory_ = "/home/hu/Desktop/EvolutionComputation/jmetalExperiment/" +
					exp.experimentName_ + "/M=" + exp.noOfObjectives_;
			exp.paretoFrontDirectory_ = "./PF";
			exp.algorithmSettings_ = new Settings[numberOfAlgorithms];

			exp.independentRuns_ = 30;

			exp.initExperiment();
//			exp.runExperiment(6);
//			exp.generateQualityIndicators();
//		exp.initExperiment();
			// Run the experiments
//		exp.runExperiment(6);
//		exp.generateQualityIndicators();
			// Generate latex tables
			// generate tables without test symbols
			exp.generateLatexTables(false);
			// generate tables with test symbols
//     exp.generateLatexTables(true) ;
			// Configure the R scripts to be generated
			int rows;
			int columns;
			String prefix;
			String[] problems;
			boolean notch;
			rows = 3;
			columns = 3;
			prefix = "DTLZ1";
			problems = new String[]{"DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7"
//				, "WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9"};
			};

//		exp.generateRBoxplotScripts(rows, columns, problems, prefix,
//				notch = false, exp);
//		exp.generateRWilcoxonScripts(problems, prefix, exp);
			// Applying Friedman test
			Friedman test = new Friedman(exp);
//    test.executeTest("EPSILON");
//		test.executeTest("HV");
//    test.executeTest("GSPREAD");
			test.executeTest("IGD");
//    test.executeTest("RUNTIME");
		}
	} // main
} // AdMOEAStudy

