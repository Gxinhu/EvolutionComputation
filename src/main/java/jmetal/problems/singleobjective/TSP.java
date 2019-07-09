//  TSP.java
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

package jmetal.problems.singleobjective;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.encodings.solutionType.PermutationSolutionType;
import jmetal.encodings.variable.Permutation;
import jmetal.util.JMException;

import java.io.*;

/**
 * Class representing a TSP (Traveling Salesman Problem) problem.
 */
public class TSP extends Problem {

	public int numberOfCities_;
	public double[][] distanceMatrix_;

	public TSP(String solutionType) {
		this(solutionType, "./TSP/eil101.tsp");
	}

	/**
	 * Creates a new TSP problem instance. It accepts data files from TSPLIB
	 *
	 * @param filename The file containing the definition of the problem
	 */
	public TSP(String solutionType, String filename) {
		numberOfVariables_ = 1;
		numberOfObjectives_ = 1;
		numberOfConstraints_ = 0;
		problemName_ = "TSP";

		length_ = new int[numberOfVariables_];
		length_[0] = numberOfCities_;

		try {
			if (solutionType.compareTo("Permutation") == 0) {
				solutionType_ = new PermutationSolutionType(this);
			} else {
				throw new JMException("Solution type invalid");
			}
		} catch (JMException e) {
			e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
		}
		try {
			readProblem(filename);
		} catch (IOException e) {
			e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
		}
		System.out.println(numberOfCities_);

	} // TSP

	/**
	 * Evaluates a solution
	 *
	 * @param solution The solution to evaluate
	 */
	public void evaluate(Solution solution) {
		double fitness;

		fitness = 0.0;

		for (int i = 0; i < (numberOfCities_ - 1); i++) {
			int x;
			int y;

			x = ((Permutation) solution.getDecisionVariables()[0]).vector_[i];
			y = ((Permutation) solution.getDecisionVariables()[0]).vector_[i + 1];
//  cout << "I : " << i << ", x = " << x << ", y = " << y << endl ;    
			fitness += distanceMatrix_[x][y];
		} // for
		int firstCity;
		int lastCity;

		firstCity = ((Permutation) solution.getDecisionVariables()[0]).vector_[0];
		lastCity = ((Permutation) solution.getDecisionVariables()[0]).vector_[numberOfCities_ - 1];
		fitness += distanceMatrix_[firstCity][lastCity];

		solution.setObjective(0, fitness);//objective value

	} // evaluate


	public void readProblem(String fileName) throws
			IOException {
		Reader inputFile = new BufferedReader(
				new InputStreamReader(
						new FileInputStream(fileName)));

		StreamTokenizer token = new StreamTokenizer(inputFile);
		try {
			boolean found;
			found = false;

			token.nextToken();
			while (!found) {
				if ((token.sval != null) && ((token.sval.compareTo("DIMENSION") == 0))) {
					found = true;
				} else {
					token.nextToken();
				}
			} // while

			token.nextToken();
			token.nextToken();

			numberOfCities_ = (int) token.nval;

			distanceMatrix_ = new double[numberOfCities_][numberOfCities_];

			//Find the EDGE_WEIGHT_TYPE
			found = false;

			token.nextToken();
			while (!found) {
				if ((token.sval != null) && ((token.sval.compareTo("TYPE") == 0))) {
					found = true;
				} else {
					token.nextToken();
				}
			} // while

			token.nextToken();
			token.nextToken();
			String type = token.sval;

			// Find the string SECTION
			found = false;
			token.nextToken();
			while (!found) {
				if ((token.sval != null) &&
						((token.sval.compareTo("SECTION") == 0))) {
					found = true;
				} else {
					token.nextToken();
				}
			} // while

			// Read the data
			if (type.compareTo("EUC") == 0) {
				double[] c = new double[2 * numberOfCities_];

				for (int i = 0; i < numberOfCities_; i++) {
					token.nextToken();
					int j = (int) token.nval;
					token.nextToken();
					double base = token.nval;

					token.nextToken();
					if ((token.sval != null) && ((token.sval.compareTo("e") == 0))) {
						token.nextToken();
						token.nextToken();
						int power = (int) token.nval;

						String str = Double.toString(base) + "e+" + Integer.toString(power);

						c[2 * (j - 1)] = Double.parseDouble(str);
//		        System.out.println(c[2*(j-1)]);        
						token.nextToken();
						token.nextToken();
						base = token.nval;

						token.nextToken();
						token.nextToken();
						power = (int) token.nval;
						str = Double.toString(base) + "e+" + Integer.toString(power);

						c[2 * (j - 1) + 1] = Double.parseDouble(str);
//		        System.out.println(c[2*(j-1)+1]);
					} else {
						c[2 * (j - 1)] = base;
//	        	  System.out.println(c[2*(j-1)]);         	
						c[2 * (j - 1) + 1] = token.nval;
//	        	 System.out.println(c[2*(j-1)+1]);
					}
				} // for

				double dist;
				for (int k = 0; k < numberOfCities_; k++) {
					distanceMatrix_[k][k] = 0;
					for (int j = k + 1; j < numberOfCities_; j++) {
						dist = Math.sqrt(Math.pow((c[k * 2] - c[j * 2]), 2.0) +
								Math.pow((c[k * 2 + 1] - c[j * 2 + 1]), 2));
						dist = (int) (dist + .5);
						distanceMatrix_[k][j] = dist;
						distanceMatrix_[j][k] = dist;
					} // for
				} // for
			} else if (type.compareTo("EXPLICIT") == 0) {
				token.nextToken();
				for (int k = 0; k < numberOfCities_; k++) {
					for (int j = 0; j < numberOfCities_; j++) {
						distanceMatrix_[k][j] = token.nval;
						System.out.print(distanceMatrix_[k][j] + " ");
						token.nextToken();
					} // for
					System.out.println("\n");
				} // for
			}

		} // try
		catch (Exception e) {
			System.err.println("TSP.readProblem(): error when reading data file " + e);
			System.exit(1);
		} // catch
	} // readProblem


	public double[][] readData(String dataType) throws JMException {
		//Always return distanceMatrix_ , since only one objective
		if (dataType.compareToIgnoreCase("distance") == 0) {
			return distanceMatrix_;
		} else if (dataType.compareToIgnoreCase("cost") == 0) {
			return distanceMatrix_;
		}
		return distanceMatrix_;
	
}

} // TSP
