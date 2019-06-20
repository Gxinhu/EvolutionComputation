/**
 * MOPSOD_main.java
 *
 * @author Xin HU
 */
package jmetal.singleobjective;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.problems.singleobjective.Rosenbrock;
import jmetal.util.JMException;

import java.util.Arrays;


public class mapsomain {
	public static void main(String[] args) throws JMException,
			SecurityException, NullPointerException, ClassNotFoundException {
		Problem problem;
		problem = new Rosenbrock("Real", 30);
//		double bias = new getcec2005bias(3).bias;
//		double bias=0;
		Algorithm algorithm;
		algorithm = new mapso(problem);
		algorithm.setInputParameter("maxIterations", 100000);
		algorithm.setInputParameter("swarmSize", 40);
		long starTime = System.currentTimeMillis();
		SolutionSet population = algorithm.execute();
		long endTime = System.currentTimeMillis();
		System.out.println(endTime - starTime + "ms");
		double minFitness = Double.POSITIVE_INFINITY;
		int minIndex = -1;
		for (int i = 0; i < population.size(); i++) {
			if (population.get(i).getObjective(0) < minFitness) {
				minFitness = population.get(i).getObjective(0);
				minIndex = i;
			}
		}
		System.out.println(population.get(minIndex).getObjective(0));
		System.out.println(Arrays.toString(population.get(minIndex).getDecisionVariables()));
	} //runtimes
} // main
