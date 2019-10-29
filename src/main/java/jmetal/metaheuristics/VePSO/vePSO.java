//  PSO.java
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

package jmetal.metaheuristics.VePSO;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.DominanceComparator;

import java.util.ArrayList;
import java.util.Comparator;

/**
 * Class implementing a single-objective PSO algorithm
 */
public class vePSO extends Algorithm {

	/**
	 * Stores the number of particles used
	 */
	private int particlesSize;
	private Comparator dominance;
	/**
	 * Stores the maximum number of iteration
	 */
	private int maxIterations;
	/**
	 * Stores the current number of iteration
	 */
	private int iteration;
	/**
	 * Stores the particles
	 */
	private SolutionSet[] particles;
	private CrowdingArchive Archive;
	/**
	 * Stores the local best solutions found so far for each particles
	 */
	private SolutionSet[] personalBest;
	/**
	 * Stores the global best solution found
	 */
	private Solution[] globalBest;
	private int numberOfSwarm;

	private int archiveSize;
	double w = 0.72;
	double c1 = 1.49;
	double c2 = 1.49;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public vePSO(Problem problem) {
		super(problem);
	} // Constructor

	/**
	 * Runs of the SMPSO algorithm.
	 *
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 * solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		maxIterations = (int) this.getInputParameter("maxIterations");
		iteration = 0;
		particlesSize = (int) this.getInputParameter("swarmSizes");
		archiveSize = (int) this.getInputParameter("archiveSize");
		Archive = new CrowdingArchive(archiveSize, problem_.getNumberOfObjectives());
		numberOfSwarm = 5;
		particles = new SolutionSet[numberOfSwarm];
		personalBest = new SolutionSet[numberOfSwarm];
		globalBest = new Solution[numberOfSwarm];
		dominance = new DominanceComparator();
		for (int i = 0; i < numberOfSwarm; i++) {
			particles[i] = new SolutionSet(particlesSize);
			personalBest[i] = new SolutionSet(particlesSize);
		}
		initSwarmAndArchive();
		while (iteration < maxIterations) {
			for (int i = 0; i < numberOfSwarm; i++) {
				updatePersonalBest(i);
				updateGlobalBest(i);
			}
			updatePopulatiton();
			iteration++;
		}
		return Archive;
	} // execute

	private void updatePersonalBest(int index) {
		for (int i = 0; i < particlesSize; i++) {
			int flag;
			flag = dominance.compare(particles[index].get(i), personalBest[index].get(i));
			if (flag == 0) {
				if (PseudoRandom.randDouble() < 0.5) {
					personalBest[index].replace(i, new Solution(particles[index].get(i)));
				}
			} else if (flag == -1) {
				personalBest[index].replace(i, new Solution(particles[index].get(i)));
			}
		}
	}

	private void updatePopulatiton() throws JMException {
		Variable[] pBest, gBest;
		double r1, r2;
		for (int swarmIndex = 0; swarmIndex < numberOfSwarm; swarmIndex++) {
			for (int i = 0; i < particlesSize; i++) {
				Variable[] particle = particles[swarmIndex].get(i).getDecisionVariables();
				double[] velocity = particles[swarmIndex].get(i).getSpeed();
				pBest = personalBest[swarmIndex].get(i).getDecisionVariables();
				int index = 0;
				if (swarmIndex == 0) {
					index = numberOfSwarm - 1;
					gBest = globalBest[index].getDecisionVariables();
				} else {
					index = swarmIndex - 1;
					gBest = globalBest[index].getDecisionVariables();
				}
				for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
					r1 = PseudoRandom.randDouble();
					r2 = PseudoRandom.randDouble();
					c1 = PseudoRandom.randDouble(1.5, 2.0);
					c2 = PseudoRandom.randDouble(1.5, 2.0);
					w = PseudoRandom.randDouble(0.1, 0.5);
					double temp = (w * velocity[j]) + c1 * r1 * (pBest[j].getValue() - particle[j].getValue())
							+ c2 * r2 * (gBest[j].getValue() - particle[j].getValue());
					particles[swarmIndex].get(i).setSpeed(j, temp);
				}
			}
			for (int n = 0; n < particlesSize; n++) {

				// DecisionVariables particle = ;
				for (int var = 0; var < particles[swarmIndex].get(n).getDecisionVariables().length; var++) {
					particles[swarmIndex].get(n).getDecisionVariables()[var].setValue(particles[swarmIndex].get(n).getDecisionVariables()[var].getValue() +
							particles[swarmIndex].get(n).getSpeed()[var]);
					if (particles[swarmIndex].get(n).getDecisionVariables()[var].getValue() < problem_.getLowerLimit(var)) {
						particles[swarmIndex].get(n).getDecisionVariables()[var].setValue(problem_.getLowerLimit(var));
					}
					if (particles[swarmIndex].get(n).getDecisionVariables()[var].getValue() > problem_.getUpperLimit(var)) {
						particles[swarmIndex].get(n).getDecisionVariables()[var].setValue(problem_.getUpperLimit(var));
					}
				}
				problem_.evaluate(particles[swarmIndex].get(n));
				Archive.add(new Solution(particles[swarmIndex].get(n)));
			}
		}
	}


	private void updateGlobalBest(int index) {
		ArrayList<Integer> swarmList = new ArrayList<>();
		for (int i = 0; i < Archive.size(); i++) {
			if (Archive.get(i).getSwarmIndex() == index) {
				swarmList.add(i);
			}
		}
		int randomIndex = 0;
		if (swarmList.size() != 0) {
			randomIndex = PseudoRandom.randInt(0, swarmList.size() - 1);
			globalBest[index] = new Solution(Archive.get(swarmList.get(randomIndex)));
		} else {
			randomIndex = PseudoRandom.randInt(0, Archive.size() - 1);
			globalBest[index] = new Solution(Archive.get(randomIndex));
		}
	}

	private void initSwarmAndArchive() throws ClassNotFoundException, JMException {
		for (int i = 0; i < numberOfSwarm; i++) {
			for (int j = 0; j < particlesSize; j++) {
				Solution solution = new Solution(problem_);
				solution.setSwarmIndex(i);
				problem_.evaluate(solution);
				particles[i].add(new Solution(solution));
				personalBest[i].add(new Solution(solution));
				Archive.add(new Solution(solution));
			}
		}
	}


} // PSO
