
package jmetal.metaheuristics.Vapso.util;

import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.util.DistanceShifted;

import java.util.Arrays;


public class ShiftedEuclideanDistanceAssigment {
	private int numberOfObjectives;

	public ShiftedEuclideanDistanceAssigment(Problem problem) {
		this.numberOfObjectives = problem.getNumberOfObjectives();

	}

	public void fitnessCompute(SolutionSet solutionSet) {
		int size = solutionSet.size();
		//Use a new SolutionSet to elite alter original solutionSet
		SolutionSet front = new SolutionSet(size);
		for (int i = 0; i < size; i++) {
			front.add(solutionSet.get(i));
		}
		normalizationByNSGA3 normalizationbynsga3 = new normalizationByNSGA3(numberOfObjectives);
		normalizationbynsga3.normalization(front);
		sdeassigntocrowdingdistancefield(front);
	}

	private void sdeassigntocrowdingdistancefield(SolutionSet solutionSet) {
		int size = solutionSet.size();
		double maxCrowdistance = -1.0e+30;
		double minCrowdistance = 1.0e+30;

		if (size == 0) {
			return;
		}
		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(1);
			return;
		} // if

		if (size == 2) {
			solutionSet.get(0).setCrowdingDistance(1);
			solutionSet.get(1).setCrowdingDistance(1);
			return;
		} // if

		for (int i = 0; i < size; i++) {
			solutionSet.get(i).setCrowdingDistance(0.0);
		}

		DistanceShifted SDEdistance_ = new DistanceShifted();
		double[][] SDEdistance = SDEdistance_.translateDistanceMatrixShifted(solutionSet);
		int t = 1;
		for (int i = 0; i < SDEdistance.length; i++) {
			Arrays.sort(SDEdistance[i]);
			double tDistance = SDEdistance[i][t];
			solutionSet.get(i).setCrowdingDistance(tDistance);
		} // for

	}


}
