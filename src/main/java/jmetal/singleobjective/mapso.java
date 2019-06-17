package jmetal.singleobjective;


import jmetal.core.*;
import jmetal.util.JMException;

public class mapso extends Algorithm {

	private Problem problem;
	private int populationSize;
	private int maxIterations;
	private Solution globalBest;
	private SolutionSet population;
	private double[][] speed;
	private SolutionSet personalBest;
	private double w, theta1, theta2;

	mapso(Problem problem) {
		super(problem);
	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		maxIterations = (Integer) this.getInputParameter("maxIterations");
		populationSize = (Integer) this.getInputParameter("swarmSize");
		population = new SolutionSet(populationSize);
		personalBest = new SolutionSet(populationSize);
		speed = new double[populationSize][problem.getNumberOfVariables()];
		initPopulation();
		globalBest = new Solution(population.get(0));
		int iteration = 0;
		while (iteration < maxIterations) {
			calculateCoefficientValues();
			updatePopulation();
			++iteration;
		}
		return population;
	}

	private void calculateCoefficientValues() {
		w = 0.5;
		theta1 = 1.5;
		theta2 = 1.5;
	}

	private void updatePopulation() throws JMException {
		Variable[] pbest, gbest;
		findGlobalbest();
		gbest = globalBest.getDecisionVariables();
		for (int i = 0; i < populationSize; i++) {
			Variable[] particle = population.get(i).getDecisionVariables();
			double[] velocity = population.get(i).getSpeed();
			if (personalBest.get(i).getFitness() > population.get(i).getFitness()) {
				personalBest.replace(i, new Solution(population.get(i)));
			}
			pbest = personalBest.get(i).getDecisionVariables();
			for (int j = 0; j < problem.getNumberOfVariables(); j++) {
				double temp = (w * velocity[j]) + theta1 * (pbest[j].getValue() - particle[j].getValue())
						+ theta2 * (gbest[j].getValue() - particle[j].getValue());
				population.get(i).setSpeed(j, temp);
			}
		}
		for (int n = 0; n < this.populationSize; n++) {

			// DecisionVariables particle = ;
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
			problem.evaluate(population.get(n));
		}
	}

	private void findGlobalbest() {
		double nowGlobalFitness = globalBest.getFitness();
		int index = -1;
		for (int i = 0; i < population.size(); i++) {
			if (population.get(i).getFitness() < nowGlobalFitness) {
				nowGlobalFitness = population.get(i).getFitness();
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
		}
	} // initPopulation

} // MOPSOD