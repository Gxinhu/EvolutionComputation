//  NSGAII.java
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

public class NSGAII extends Algorithm {

	private double[][] realtimeIGD;
	private double[][] realtimeSpeard;
	private double[][] realtimeGD;
	int evaluations;
	QualityIndicator indicators;
	boolean save;
	int runtimes;
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public NSGAII(Problem problem, boolean save, int runtimes) {
		super(problem);
		this.save = save;
		this.runtimes = runtimes;
	} // NSGAII

	/**
	 * Runs the NSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int populationSize;
		int maxEvaluations;
		// QualityIndicator object
		int requiredEvaluations; // Use in the example of use of the
		// indicators object (see below)

		SolutionSet population;
		SolutionSet offspringPopulation;
		SolutionSet union;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		Distance distance = new Distance();

		// Read the parameters
		populationSize = (Integer) getInputParameter("swarmSize");
		maxEvaluations = (Integer) getInputParameter("maxIterations");
		indicators = (QualityIndicator) getInputParameter("indicators");
		realtimeIGD = new double[maxEvaluations / 10 + 1][2];
		realtimeSpeard = new double[maxEvaluations / 10 + 1][2];
		realtimeGD = new double[maxEvaluations / 10 + 1][2];
		// Initialize the variables
		population = new SolutionSet(populationSize);
		evaluations = 0;

		requiredEvaluations = 0;

		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");

		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			//interation++;
			population.add(newSolution);
		} // for
		Ranking rankings = new Ranking(population);
		calulateindicator(rankings.getSubfront(0));
		// population.printFeasibleFUN("initialsb_NSGAII");
		// Generations
		while (evaluations < maxEvaluations) {

			// Create the offSpring solutionSet
			offspringPopulation = new SolutionSet(populationSize);
			Solution[] parents = new Solution[2];
			for (int i = 0; i < (populationSize / 2); i++) {
				if (evaluations < maxEvaluations) {
					// obtain parents
					parents[0] = (Solution) selectionOperator
							.execute(population);
					parents[1] = (Solution) selectionOperator
							.execute(population);
					Solution[] offSpring = (Solution[]) crossoverOperator
							.execute(parents);
					mutationOperator.execute(offSpring[0]);
					mutationOperator.execute(offSpring[1]);
					problem_.evaluate(offSpring[0]);
					problem_.evaluateConstraints(offSpring[0]);
					problem_.evaluate(offSpring[1]);
					problem_.evaluateConstraints(offSpring[1]);
					offspringPopulation.add(offSpring[0]);
					offspringPopulation.add(offSpring[1]);
					//interation += 2;
				} // if
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


			evaluations++;
			if (evaluations % 10 == 0) {
				rankings = new Ranking(population);
				calulateindicator(rankings.getSubfront(0));
			}

		} // while


		// Return as output parameter the required interation
		//setOutputParameter("interation", requiredEvaluations);
		if (this.save) {
			savetofile savetofile = new savetofile(problem_, "out/nsga2/IGD/" + problem_.getName(), runtimes, realtimeIGD);
			savetofile.save();
			savetofile = new savetofile(problem_, "out/nsga2/GD/" + problem_.getName(), runtimes, realtimeGD);
			savetofile.save();
			savetofile = new savetofile(problem_, "out/nsga2/Spread/" + problem_.getName(), runtimes, realtimeSpeard);
			savetofile.save();
		}
		// Return the first non-dominated front
		Ranking ranking = new Ranking(population);
//		ranking.getSubfront(0).printFeasibleFUN("FUNsb_NSGAII");

		return ranking.getSubfront(0);
	} // execute
	private void calulateindicator(SolutionSet archive) {
		if (this.save) {
			realtimeSpeard[evaluations / 10][0] = evaluations;
			realtimeSpeard[evaluations / 10][1] = indicators.getGeneralizedSpread(archive);
			realtimeIGD[evaluations / 10][0] = evaluations;
			realtimeIGD[evaluations / 10][1] = indicators.getCEC_IGD(archive);
			realtimeGD[evaluations / 10][0] = evaluations;
			realtimeGD[evaluations / 10][1] = indicators.getGD(archive);
		}
	}
} // NSGA-II
