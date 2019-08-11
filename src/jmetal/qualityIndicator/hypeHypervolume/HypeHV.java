package jmetal.qualityIndicator.hypeHypervolume;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;

import java.io.IOException;

/**
 * this is calualate HV in RVEA
 */
public class HypeHV {
	public jmetal.qualityIndicator.util.MetricsUtil utils_;
	private double[][] pf_;
	private double[][] pfMatrix_ = null;
	private String problemNames = null;

	public HypeHV(double[][] paretoFront, double[][] pfMatrix) {
		pf_ = paretoFront;
		pfMatrix_ = pfMatrix;
		utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	} // Constructor

	public HypeHV(double[][] solutionFront, String problemName) {
		problemNames = problemName;
		pf_ = solutionFront;
		utils_ = new jmetal.qualityIndicator.util.MetricsUtil();
	}


	public static double hv2point(Solution Point1, Solution ref) {
		double x = ref.getObjective(0) - Point1.getObjective(0);
		for (int j = 1; j < Point1.getNumberOfObjectives(); j++) {
			x = x * (ref.getObjective(j) - Point1.getObjective(j));
		}
		return x;
	}

	public double calculatewfghv() throws IOException {
		SolutionSet sb = new SolutionSet(pf_.length);
		for (int ss = 0; ss < pf_.length; ss++) {
			Solution sss = new Solution(pf_[ss].length);
			for (int j = 0; j < pf_[ss].length; j++) {

				sss.setObjective(j, pf_[ss][j]);
			}
			sb.add(sss);


		}

		double hv = 0;
		int number = pf_[0].length;
		Solution referencePoint1 = new Solution(number);
		double[] maxObjectives = new double[number];
		for (int i = 0; i < number; i++) {
			maxObjectives[i] = 0;
		}
		if (problemNames != null) {
			for (int j = 0; j < number; j++) {
				if ("DTLZ1".equals(problemNames)) {
					referencePoint1.setObjective(j, 0.5);
				} else if (problemNames.equals("DTLZ2") | problemNames.equals("DTLZ3") | problemNames.equals("DTLZ4")) {
					referencePoint1.setObjective(j, 1.0);
				} else if (problemNames.equals("DTLZ5") | problemNames.equals("DTLZ6")) {
					if (j != number - 1) {
						referencePoint1.setObjective(number - j - 1, Math.pow(Math.sqrt(2) / 2, j));
					} else {
						referencePoint1.setObjective(0, referencePoint1.getObjective(1));
					}
				} else if (problemNames.equals("DTLZ7")) {
					if (j != number - 1) {
						referencePoint1.setObjective(j, 0.8594);
					} else {
						referencePoint1.setObjective(j, 2 * number);
					}
				} else if (problemNames.contains("WFG")) {
					referencePoint1.setObjective(j, 2.0 * (j + 1));
				}
			}

		} else {
			for (double[] doubles : pfMatrix_) {
				for (int j = 0; j < number; j++) {
					if (maxObjectives[j] < doubles[j]) {
						maxObjectives[j] = doubles[j];
					}
				}
			}
			for (int i = 0; i < number; i++) {
				referencePoint1.setObjective(i, maxObjectives[i]);
			}
		}
		//NORMALIZATION
		double[] minimumValues = utils_.getMinimumValues(pf_, number);
		for (int i = 0; i < number; i++) {
			if (minimumValues[i] >= 0) {
				minimumValues[i] = 0;
			}
		}
		for (int j = 0; j < pf_.length; j++) {
			for (int k = 0; k < number; k++) {
				sb.get(j).setObjective(k, (sb.get(j).getObjective(k) - minimumValues[k]) / (1.1 * referencePoint1.getObjective(k) - minimumValues[k]));
			}
		}
		SolutionSet invertedFront;
		//invertedFront = utils_.invertedFront(pf_,number);
		for (int j = 0; j < pf_.length; j++) {
			for (int k = 0; k < number; k++) {
				if (sb.get(j).getObjective(k) - 1e-10 > (1.0)) {
					//pf_.remove(j);
					//j--;
					for (int s = 0; s < number; s++) {
						sb.get(j).setObjective(s, 1.0);
					}
					break;
				}
			}
		}
		if (sb.size() == 0) {
			return 0;
		}
		for (int j = 0; j < number; j++)
		//referencePoint1.setObjective(j,2.0*j+10);
		{
			referencePoint1.setObjective(j, 1.0);
		}
		if (sb.size() == 0) {
			hv = 0.0;
		} else if (sb.size() == 1) {
			hv = hv2point(sb.get(0), referencePoint1);
		} else if (sb.size() == 2) {
			double[] mid = new double[number];
			for (int j = 0; j < number; j++) {
				mid[j] = Math.max(sb.get(0).getObjective(j), sb.get(1).getObjective(j));
			}
			Solution midp = new Solution(number);
			for (int i = 0; i < number; i++) {
				midp.setObjective(i, mid[i]);
			}
			hv = hv2point(sb.get(0), referencePoint1) + hv2point(sb.get(1), referencePoint1) - hv2point(midp, referencePoint1);

		} else if (sb.size() == 3) {
			double[] w01 = new double[number];
			double[] w02 = new double[number];
			double[] w12 = new double[number];
			double[] w012 = new double[number];
			for (int j = 0; j < number; j++) {
				w01[j] = Math.max(sb.get(0).getObjective(j), sb.get(1).getObjective(j));
				w02[j] = Math.max(sb.get(0).getObjective(j), sb.get(2).getObjective(j));
				w12[j] = Math.max(sb.get(1).getObjective(j), sb.get(2).getObjective(j));
			}
			for (int j = 0; j < number; j++) {
				w012[j] = Math.max(w02[j], sb.get(1).getObjective(j));
			}
			Solution p01 = new Solution(number);
			Solution p02 = new Solution(number);
			Solution p12 = new Solution(number);
			Solution p012 = new Solution(number);
			for (int i = 0; i < number; i++) {
				p01.setObjective(i, w01[i]);
				p02.setObjective(i, w02[i]);
				p12.setObjective(i, w12[i]);
				p012.setObjective(i, w012[i]);
			}
			hv = hv2point(sb.get(0), referencePoint1) + hv2point(sb.get(1), referencePoint1) + hv2point(sb.get(2), referencePoint1)
					- hv2point(p01, referencePoint1) - hv2point(p02, referencePoint1) - hv2point(p12, referencePoint1) + hv2point(p012, referencePoint1);
		} else {
			HypEHypervolume hype = new HypEHypervolume();
			hv = hype.estimateHypervolume(sb, referencePoint1, 1000000);
		}
		return hv;
	} // CalculateHypervolume
}
