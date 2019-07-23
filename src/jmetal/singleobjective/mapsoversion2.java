package jmetal.singleobjective;


import jmetal.core.*;
import jmetal.util.JMException;

import java.util.Random;

public class mapsoversion2 extends Algorithm {

	private Problem problem;
	private int populationSize;
	private int maxIterations;
	private Solution globalBest;
	private SolutionSet population, lastPopulation;
	private double[][] speed;
	private Random random = new Random();
	private SolutionSet personalBest;
	private double w, theta1, theta2;
	private int iteration1, iteration2, iteration;
	private double vMax, vMin, pMin, pMax, fMax, fMin, mu1, mu2, sigma1, sigma2;

	mapsoversion2(Problem problem) {
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
		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		population = new SolutionSet(populationSize);
		lastPopulation = new SolutionSet(populationSize);
		personalBest = new SolutionSet(populationSize);
		iteration1 = maxIterations / 5;
		iteration2 = 4 * maxIterations / 5;
		initPopulation();
		globalBest = new Solution(population.get(0));
		findGlobalbest();
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
		} else if (iteration > iteration1 & iteration < (iteration2 - iteration1) / 2) {
			p1 = (iteration - iteration1) * (pMax - pMin) / ((double) (iteration2 - iteration1) / 2 - iteration1) + pMin;
		} else if (iteration < iteration2 & iteration > (iteration2 - iteration1) / 2) {
			p1 = (iteration - iteration1) * (pMin - pMax) / (iteration2 - (double) (iteration2 - iteration1) / 2) + pMax;
		} else {
			p1 = pMin;
		}
		//calculate F
		double F;
		if (iteration < iteration1) {
			F = fMin;
		} else if (iteration > iteration1 & iteration < iteration2) {
			F = 1;
		} else {
			F = fMax;
		}
		double a = Math.sqrt(F);
		double m1 = (a + 1) * (a + 1) * (a * a + 3 * a + 1);
		double m2 = (a + 1) * (a + 1) * (2 * a * a + 3 * a + 2);
		w = (m1 * vc + m2 * p1 * vc + p1 - 1) / (m2 * vc + m1 * p1 * vc - p1 + 1);
		double c = 2 * (1 - p1) * (w + 1) / (a + 1);
		mu1 = c / 2;
		mu2 = a * c / 2;
		sigma1 = c / Math.sqrt(12);
		sigma2 = a * c / Math.sqrt(12);
	}

	private void updatePopulation() throws JMException {
		Variable[] pbest, gbest;
		gbest = globalBest.getDecisionVariables();

		for (int i = 0; i < populationSize; i++) {
			theta1 = random.nextGaussian() * sigma1 + mu1;
			theta2 = random.nextGaussian() * sigma2 + mu2;
			double l = 1 + w - theta1 - theta2;
			Variable[] particle = population.get(i).getDecisionVariables();
			Variable[] lastParticle = lastPopulation.get(i).getDecisionVariables();
			lastPopulation.replace(i, new Solution(population.get(i)));

			pbest = personalBest.get(i).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {

				double temp = l * particle[j].getValue() - w * lastParticle[j].getValue() + theta1 * pbest[j].getValue()
						+ theta2 * gbest[j].getValue();
				population.get(i).getDecisionVariables()[j].setValue(temp);
			}

		}
		for (int n = 0; n < this.populationSize; n++) {
			// DecisionVariables particle ;
			for (int var = 0; var < population.get(n).getDecisionVariables().length; var++) {
				if (population.get(n).getDecisionVariables()[var].getValue() < problem.getLowerLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getLowerLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
				if (population.get(n).getDecisionVariables()[var].getValue() > problem.getUpperLimit(var)) {
					population.get(n).getDecisionVariables()[var].setValue(problem.getUpperLimit(var));
					population.get(n).setSpeed(var, 0.0);
				}
			}
			problem.evaluate(population.get(n));
			if (personalBest.get(n).getObjective(0) > population.get(n).getObjective(0)) {
				personalBest.replace(n, new Solution(population.get(n)));
			}
		}
		findGlobalbest();
	}

	private void findGlobalbest() {
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
			lastPopulation.add(new Solution(newSolution));
		}
	} // initPopulation

} // MOPSOD