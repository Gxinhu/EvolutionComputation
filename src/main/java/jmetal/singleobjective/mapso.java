/*
  mapso.java

  @author Xin HU
 * A theoretical guideline for designing an effective adaptive particle swarm.
 */
package jmetal.singleobjective;
import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.parallel.MultithreadedEvaluator;

import java.util.Random;

public class mapso extends Algorithm {

	private Problem problem;
	private int populationSize;
	private Solution globalBest;
	private Random random = new Random();
	private SolutionSet population;
	private SolutionSet personalBest;
	private double w;
	private int iteration1, iteration2, iteration;
	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;

	mapso(Problem problem) {
		super(problem);
		this.problem = problem;
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		//init coefficient
		vMax = 25.0;
		vMin = 5.0;
		pMax = 0.8;
		pMin = 0.1;
		fMax = 25.0;
		fMin = 0.25;
		int maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		population = new SolutionSet(populationSize);
		personalBest = new SolutionSet(populationSize);
		iteration1 = maxIterations / 5;
		iteration2 = 4 * maxIterations / 5;
		initPopulation();
		globalBest = new Solution(population.get(0));
		findGlobalBest();
		//Start
		iteration = 0;
		while (iteration < maxIterations) {
			calculateCoefficientValues();
			updatePopulation();
			++iteration;
//
		}
		return population;
	}

	private void calculateCoefficientValues() {
		//calculate Vc
		double vc;
		if (iteration < iteration1) {
			vc = vMax;
		} else if (iteration > iteration1 & iteration < iteration2) {
			vc = (iteration - iteration1) * (vMin - vMax) / (iteration2 - iteration1) + vMax;
		} else {
			vc = vMin;
		}
		//calculate p1
		double p1;
		if (iteration < iteration1) {
			p1 = pMin;
		} else if (iteration > iteration1 & iteration < (iteration2 - iteration1) / 2.0) {
			p1 = (iteration - iteration1) * (pMax - pMin) / ((double) (iteration2 - iteration1) / 2.0 - iteration1) + pMin;
		} else if (iteration < iteration2 & iteration > (iteration2 - iteration1) / 2.0) {
			p1 = (iteration - iteration1) * (pMin - pMax) / (iteration2 - (double) (iteration2 - iteration1) / 2.0) + pMax;
		} else {
			p1 = pMin;
		}
		//calculate F
		double f;
		if (iteration < iteration1) {
			f = fMin;
		} else if (iteration > iteration1 & iteration < iteration2) {
			f = 1;
		} else {
			f = fMax;
		}
		double a = Math.sqrt(f);
		double m1 = (a + 1.0) * (a + 1.0) * (a * a + 3.0 * a + 1.0);
		double m2 = (a + 1.0) * (a + 1.0) * (2.0 * a * a + 3.0 * a + 2.0);
		w = (m1 * vc + m2 * p1 * vc + p1 - 1.0) / (m2 * vc + m1 * p1 * vc - p1 + 1.0);
		double c = 2.0 * (1 - p1) * (w + 1.0) / (a + 1.0);
		mu1 = c / 2.0;
		mu2 = a * c / 2.0;
		sigma1 = c / Math.sqrt(12.0);
		sigma2 = a * c / Math.sqrt(12.0);
	}

	private void updatePopulation() throws JMException {
		Variable[] pbest, gbest;
		gbest = globalBest.getDecisionVariables();
		for (int i = 0; i < populationSize; i++) {
//			theta1 = random.nextGaussian() * sigma1 + mu1;
//			theta2 = random.nextGaussian() * sigma2 + mu2;
			Variable[] particle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			pbest = personalBest.get(i).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				double theta1 = random.nextGaussian() * sigma1 + mu1;
				double theta2 = random.nextGaussian() * sigma2 + mu2;
				double temp = (w * velocity[j]) + theta1 * (pbest[j].getValue() - particle[j].getValue())
						+ theta2 * (gbest[j].getValue() - particle[j].getValue());
				population.get(i).setSpeed(j, temp);
			}
		}
		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle ;
			for (int var = 0; var < population.get(n).getDecisionVariables().length; var++) {
				population.get(n).getDecisionVariables()[var].setValue(population.get(n).getDecisionVariables()[var].getValue() + population.get(n).getSpeed()[var]);
				if (population.get(n).getDecisionVariables()[var].getValue() < problem.getLowerLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getLowerLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
				if (population.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
			}
//			problem.evaluate(population.get(n));

		}
		MultithreadedEvaluator multithreadedEvaluator = new MultithreadedEvaluator(populationSize);
		multithreadedEvaluator.startEvaluator(problem);
		for (int n = 0; n < populationSize; n++) {
			multithreadedEvaluator.addSolutionForEvaluation(population.get(n));
		}
		multithreadedEvaluator.parallelEvaluation();
		multithreadedEvaluator.stopEvaluator();
		findGlobalBest();
	}

	private void findGlobalBest() {
		double nowGlobalFitness = globalBest.getObjective(0);
		int index = -1;
		for (int i = 0; i < population.size(); i++) {
			if (population.get(i).getObjective(0) < nowGlobalFitness) {
				nowGlobalFitness = population.get(i).getObjective(0);
				index = i;
			}
		}
		if (index != -1) {
			globalBest = new Solution(population.get(index));
		}
	}


	public void initPopulation() throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem);
			problem.evaluate(newSolution);
			population.add(newSolution);
			personalBest.add(new Solution(newSolution));
		}
	} // initPopulation

} // MOPSOD