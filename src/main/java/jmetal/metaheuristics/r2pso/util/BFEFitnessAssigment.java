//  BFEFitnessAssigment.java
//
//  Author:
//      S.B. Liu <2150230420@email.szu.edu.cn>
//  Copyright (c) 2016 S.B. Liu 
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
package jmetal.metaheuristics.r2pso.util;

import jmetal.core.SolutionSet;
import jmetal.util.DistanceShifted;
import jmetal.util.comparators.ObjectiveComparator;

import java.util.Arrays;

/**
 * This class defined a balanceable fitness estimation method(as
 * used in NMPSO).
 */
public class BFEFitnessAssigment {
	private int numberOfObjectives;
	private double a, b;
	private double[] d1, d2;
	private double aveD1, aveD2;
	private double aveCrowdingDistance;
	private double[] objectiveValue;
	private double aveObjectiveValue;

	public BFEFitnessAssigment(int numberOfObjectives) {
		this.numberOfObjectives = numberOfObjectives;

		aveCrowdingDistance = 0.0;
		aveObjectiveValue = 0.0;
	}


	private void SDEAssignToCrowdingDistanceField(SolutionSet solutionSet) {
		int size = solutionSet.size();
		double maxCrowdistance = -1.0e+30;
		double minCrowdistance = 1.0e+30;

		if (size == 0) {
			return;
		}
		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
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

		for (int n = 0; n < size; n++) {
			if (solutionSet.get(n).getCrowdingDistance() > maxCrowdistance) {
				maxCrowdistance = solutionSet.get(n).getCrowdingDistance();
			}
			if (solutionSet.get(n).getCrowdingDistance() < minCrowdistance) {
				minCrowdistance = solutionSet.get(n).getCrowdingDistance();
			}
		}

		//normalize
		double sumCrowdingDistance = 0.0;
		aveCrowdingDistance = 0.0;


		for (int n = 0; n < size; n++) {
			solutionSet.get(n).setCrowdingDistance((solutionSet.get(n).getCrowdingDistance() - minCrowdistance) / (maxCrowdistance - minCrowdistance));
			sumCrowdingDistance += solutionSet.get(n).getCrowdingDistance();
		}

		aveCrowdingDistance = sumCrowdingDistance / size;
	}

	private void normalizePopulation(SolutionSet solutionSet) {
		// Obtains the lower and upper bounds of the population
		double[] maximumValues = new double[numberOfObjectives];
		double[] minimumValues = new double[numberOfObjectives];

		for (int i = 0; i < numberOfObjectives; i++) {
			maximumValues[i] = -Double.MAX_VALUE; // i.e., the minus maxium value
			minimumValues[i] = Double.MAX_VALUE; // i.e., the maximum value
		}
	   /* for (int pos = 0; pos < solutionSet.size(); pos++) {
	        for (int obj = 0; obj < numberOfObjectives; obj++) {
	          double value = solutionSet.get(pos).getObjective(obj);
	          if (value > maximumValues[obj])
	              maximumValues[obj] = value;
	          if (value < minimumValues[obj])
	              minimumValues[obj] = value;
	        }
	    }*/

		for (int i = 0; i < numberOfObjectives; i++) {
			// Sort the population by Obj n
			solutionSet.sort(new ObjectiveComparator(i));
			minimumValues[i] = solutionSet.get(0).getObjective(i);   //minimumValues
			maximumValues[i] = solutionSet.get(solutionSet.size() - 1).getObjective(i); //maximumValues
		}
		/**
		 * normalize the solutions objectives
		 * */
		for (int pos = 0; pos < solutionSet.size(); pos++) {
			for (int obj = 0; obj < numberOfObjectives; obj++) {
				double val = (solutionSet.get(pos).getObjective(obj) - minimumValues[obj]) / (maximumValues[obj] - minimumValues[obj]);
				if (val < 0 || val > 1) {
					System.out.println("val == " + val);
				}
				solutionSet.get(pos).setTranslatedObjectives(val, obj);
			}
		}

	}

}
