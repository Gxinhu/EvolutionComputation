package jmetal.metaheuristics.nsgaII;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.savesomething.savetofile;

/**
 * Implementation of NSGA-II. This implementation of NSGA-II makes use of a
 * QualityIndicator object to obtained the convergence speed of the algorithm.
 * This version is used in the paper: A.J. Nebro, J.J. Durillo, C.A. Coello
 * Coello, F. Luna, E. Alba
 * "A Study of Convergence Speed in Multi-Objective Metaheuristics." To be
 * presented in: PPSN'08. Dortmund. September 2008.
 */

public class NSGA_DE extends Algorithm {
	private double[][] realtimeIGD;
	private double[][] realtimeSpeard;
	private double[][] realtimeGD;
	int interation;
	QualityIndicator indicators;
	boolean save;
	int runtimes;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public NSGA_DE(Problem problem, boolean save, int runtimes) {
		super(problem);
		this.save = save;
		this.runtimes = runtimes;
	}

	/**
	 * Runs the NSGA-II algorithm.
	 *
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 * solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int populationSize;
		int maxIntersion;
		int requiredEvaluations; // Use in the example of use of the
		// indicators object (see below)
		int neighborsize = 20;
		int[][] neighborhood;

		SolutionSet population;
		SolutionSet offspringPopulation;
		SolutionSet union;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		Distance distance = new Distance();
		neighborhood = new int[problem_.getNumberOfObjectives()][neighborsize];

		// Read the parameters
		populationSize = (Integer) getInputParameter("swarmSize");
		maxIntersion = (Integer) getInputParameter("maxIterations");
		indicators = (QualityIndicator) getInputParameter("indicators");
		realtimeIGD = new double[maxIntersion / 10 + 1][2];
		realtimeSpeard = new double[maxIntersion / 10 + 1][2];
		realtimeGD = new double[maxIntersion / 10 + 1][2];
		// Initialize the variables
		population = new SolutionSet(populationSize);
		interation = 0;

		requiredEvaluations = 0;

		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");
		Distance distance_ = new Distance();

		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population.add(newSolution);
		} // for
		Ranking rankings = new Ranking(population);
		calulateindicator(rankings.getSubfront(0));
		interation++;
		// Generations
		while (interation < maxIntersion) {

			// Create the offSpring solutionSet
			offspringPopulation = new SolutionSet(populationSize);
			Solution[] parents = new Solution[4];
			for (int i = 0; i < populationSize; i++) {
				// obtain parents
				parents[2] = (Solution) selectionOperator
						.execute(population);
				parents[0] = (Solution) selectionOperator
						.execute(population);
				parents[1] = (Solution) selectionOperator
						.execute(population);
				double mindistance = distance_.distanceBetweenSolutions(parents[0], parents[2]);
				double mindistance2 = distance_.distanceBetweenSolutions(parents[1], parents[2]);
				if (mindistance > mindistance2) {
					parents[0] = parents[1];
				}
				parents[3] = (Solution) selectionOperator
						.execute(population);
				parents[1] = (Solution) selectionOperator
						.execute(population);
				mindistance = distance_.distanceBetweenSolutions(parents[3], parents[2]);
				mindistance2 = distance_.distanceBetweenSolutions(parents[1], parents[2]);
				if (mindistance < mindistance2) {
					parents[1] = parents[3];
				}
				Solution child = (Solution) crossoverOperator.execute(new Object[]{
						parents[2], parents});
				mutationOperator.execute(child);
				//mutationOperator.execute(offSpring[1]);
				problem_.evaluate(child);
				offspringPopulation.add(child);
			} // for

			// Create the solutionSet union of solutionSet and offSpring
			union = ((SolutionSet) population).union(offspringPopulation);

			// Ranking the union
			Ranking ranking = new Ranking(union);

			int remain = populationSize;
			int index = 0;
			SolutionSet front = null;
			population.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			while ((remain > 0) && (remain >= front.size())) {
				// Assign crowding distance to individuals
				distance.crowdingDistanceAssignment(front,
						problem_.getNumberOfObjectives());
				// Add the individuals of this front
				for (int k = 0; k < front.size(); k++) {
					population.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			} // while

			// Remain is less than front(index).size, insert only the best one
			if (remain > 0) { // front contains individuals to insert
				distance.crowdingDistanceAssignment(front,
						problem_.getNumberOfObjectives());
				front.sort(new CrowdingComparator());
				for (int k = 0; k < remain; k++) {
					population.add(front.get(k));
				} // for

				remain = 0;
			} // if
			interation++;
			if (interation % 10 == 0) {
				rankings = new Ranking(population);
				calulateindicator(rankings.getSubfront(0));
			}


		} // while

		// Return as output parameter the required interation
//		setOutputParameter("interation", requiredEvaluations);
		if (this.save) {
			savetofile savetofile = new savetofile(problem_, "out/nsga2de/IGD/" + problem_.getName(), runtimes, realtimeIGD);
			savetofile.save();
			savetofile = new savetofile(problem_, "out/nsga2de/GD/" + problem_.getName(), runtimes, realtimeGD);
			savetofile.save();
			savetofile = new savetofile(problem_, "out/nsga2de/Spread/" + problem_.getName(), runtimes, realtimeSpeard);
			savetofile.save();
		}
		// Return the first non-dominated front
		Ranking ranking = new Ranking(population);
		return ranking.getSubfront(0);
	} // execute

	private void calulateindicator(SolutionSet archive) {
		if (this.save) {
			realtimeSpeard[interation / 10][0] = interation;
			realtimeSpeard[interation / 10][1] = indicators.getGeneralizedSpread(archive);
			realtimeIGD[interation / 10][0] = interation;
			realtimeIGD[interation / 10][1] = indicators.getCEC_IGD(archive);
			realtimeGD[interation / 10][0] = interation;
			realtimeGD[interation / 10][1] = indicators.getGD(archive);
		}
	}
} // NSGA-II

