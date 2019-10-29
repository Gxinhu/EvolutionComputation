//  Rosenbrock.java
//
//  Author:
//       Esteban López-Camacho <esteban@lcc.uma.es>
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

package jmetal.problems.singleobjective.cec2005.funcition;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.problems.singleobjective.cec2005.benchmark;
import jmetal.problems.singleobjective.cec2005.test_func;
import jmetal.util.JMException;

public class cec2005F19 extends Problem {

	/**
	 * Constructor
	 * Creates a default instance of the Rosenbrock problem
	 *
	 * @param numberOfVariables Number of variables of the problem
	 * @param solutionType      The solution type must "Real" or "BinaryReal".
	 */
	public cec2005F19(String solutionType, Integer numberOfVariables) throws ClassNotFoundException {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = 1;
		numberOfConstraints_ = 0;
		problemName_ = "F19_rotated_hybrid_composition_2_narrow_basin_global_opt";

		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < numberOfVariables_; var++) {
			lowerLimit_[var] = -5;
			upperLimit_[var] = 5;
		} // for

		if (solutionType.compareTo("BinaryReal") == 0) {
			solutionType_ = new BinaryRealSolutionType(this);
		} else if (solutionType.compareTo("Real") == 0) {
			solutionType_ = new RealSolutionType(this);
		} else {
			System.out.println("Error: solution type " + solutionType + " invalid");
			System.exit(-1);
		}

	} // Rosenbrock

	/**
	 * Evaluates a solution
	 *
	 * @param solution The solution to evaluate
	 * @throws JMException
	 */
	@Override
	public void evaluate(Solution solution) throws JMException {
		Variable[] decisionVariables = solution.getDecisionVariables();

		double sum = 0.0;
		double[] x = new double[numberOfVariables_];

		for (int i = 0; i < numberOfVariables_; i++) {
			x[i] = decisionVariables[i].getValue();
		}
		// Create a benchmark object
		benchmark theBenchmark = new benchmark();
		// Use the factory function call to create a test function object
		//		test function 3 with 50 dimension
		//		the object class is "test_func"
		test_func aTestFunc = theBenchmark.testFunctionFactory(19, numberOfVariables_);
		// Invoke the function with x
		double result = aTestFunc.f(x);

		solution.setObjective(0, result);
	} // evaluate

} // Rosenbrock

