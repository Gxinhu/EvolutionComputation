//  SPEA2.java
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

package jmetal.metaheuristics.SPEA2_SDE;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Ranking;
import jmetal.util.Spea2_SDEFitness;

/**
 * This class representing the SPEA2 algorithm
 */
public class SPEA2_SDE_DE extends Algorithm {

	/**
	 * Defines the number of tournaments for creating the mating pool
	 */
	public static final int TOURNAMENTS_ROUNDS = 1;

	/**
	 * Constructor.
	 * Create a new SPEA2 instance
	 *
	 * @param problem Problem to solve
	 */
	public SPEA2_SDE_DE(Problem problem) {
		super(problem);
	} // Spea2

	/**
	 * Runs of the Spea2 algorithm.
	 *
	 * @return a <code>SolutionSet</code> that is a set of non dominated solutions
	 * as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int populationSize, archiveSize, maxEvaluations, evaluations;
		Operator crossoverDE_, mutationOperator, selectionOperator;
		SolutionSet solutionSet, archive, offSpringSolutionSet;

		//Read the params
		populationSize = ((Integer) getInputParameter("populationSize")).intValue();
		archiveSize = ((Integer) getInputParameter("archiveSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

		//Read the operators
		crossoverDE_ = operators_.get("crossoverDE"); // set the crossover operator
		mutationOperator = operators_.get("mutation");
		selectionOperator = operators_.get("selection");

		//Initialize the variables
		solutionSet = new SolutionSet(populationSize);
		archive = new SolutionSet(archiveSize);
		evaluations = 0;

		//-> Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			evaluations++;
			solutionSet.add(newSolution);
		}

		while (evaluations < maxEvaluations) {
			SolutionSet union = solutionSet.union(archive);
			Spea2_SDEFitness spea = new Spea2_SDEFitness(union);
			spea.fitnessAssign();
			archive = spea.environmentalSelection(archiveSize);
			// Create a new offspringPopulation
			offSpringSolutionSet = new SolutionSet(populationSize);

			for (int i = 0; i < archive.size(); i++) {
				Solution child;
				int r1, r2, r3;
				do {
					r1 = PseudoRandom.randInt(0, archive.size() - 1);
				} while (r1 == i);
				do {
					r2 = PseudoRandom.randInt(0, archive.size() - 1);
				} while (r2 == i || r2 == r1);
				do {
					r3 = PseudoRandom.randInt(0, archive.size() - 1);
				} while (r3 == i || r3 == r2 || r3 == r1);
				Solution[] parents = new Solution[3];
				parents[0] = archive.get(r1);
				parents[1] = archive.get(r2);
				parents[2] = archive.get(r3);
				child = (Solution) crossoverDE_.execute(new Object[]{archive.get(i), parents});
				mutationOperator.execute(child);
				problem_.evaluate(child);
				problem_.evaluateConstraints(child);
				offSpringSolutionSet.add(child);
				evaluations++;
			} // for
			// End Create a offSpring solutionSet
			solutionSet = offSpringSolutionSet;
		} // while

		Ranking ranking = new Ranking(archive);
		return ranking.getSubfront(0);
	} // execute
} // SPEA2
