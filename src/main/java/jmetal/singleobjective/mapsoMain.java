/**
 * MOPSOD_main.java
 *
 * @author Noura Al Moubayed
 */
package jmetal.singleobjective;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.problems.singleObjective.Rosenbrock;
import jmetal.util.JMException;

import java.io.IOException;


public class mapsoMain {
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		Problem problem;
		problem = new Rosenbrock("Real", 30);
		Algorithm algorithm;
		algorithm = new mapso(problem);
		algorithm.setInputParameter("maxIterations", 2000);
		algorithm.setInputParameter("swarmSize", 100);
		long startime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		long endtime = System.currentTimeMillis();

	} //runtimes
} // main
