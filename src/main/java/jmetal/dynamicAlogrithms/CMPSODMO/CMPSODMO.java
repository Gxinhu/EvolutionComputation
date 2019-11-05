//  CMPSODMO.java
//
//  Author: Xin Hu <laohuxin@gmail.com>
// Reference  [1] R. Liu, J. Li, J. fan, C. Mu, and L. Jiao, “A coevolutionary technique based on multi-swarm particle swarm optimization for dynamic multi-objective optimization,”
// Eur. J. Oper. Res., vol. 261, no. 3, pp. 1028–1051, 2017.
//


package jmetal.dynamicAlogrithms.CMPSODMO;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.savesomething.savetofile;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Comparator;

/**
 * Class implementing a single-objective PSO algorithm
 */
public class CMPSODMO extends Algorithm {

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
	private savetofile save;
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
	private SolutionSet globalBest;
	private int numberOfSwarm;

	private int archiveSize;
	double wFinal = 0.4;
	double wInitial = 0.9;
	double c1, c2, c3 = 4.0 / 3.0;
	double threshold = 0.001;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public CMPSODMO(Problem problem) {
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
		numberOfSwarm = problem_.getNumberOfObjectives();
		particles = new SolutionSet[numberOfSwarm];
		personalBest = new SolutionSet[numberOfSwarm];
		globalBest = new SolutionSet(numberOfSwarm);
		dominance = new DominanceComparator();
		for (int i = 0; i < numberOfSwarm; i++) {
			particles[i] = new SolutionSet(particlesSize);
			personalBest[i] = new SolutionSet(particlesSize);
			globalBest.add(new Solution());
		}
		initSwarmAndArchive();
		for (int i = 0; i < numberOfSwarm; i++) {
			updateGlobalBest(i);
		}
		while (iteration < maxIterations) {
			if (iteration % 10 == 0) {
				save = new savetofile(problem_, "./draw" + "/CMPSODMO" + "/PF", iteration, problem_.getPF());
				save.save();
				Archive.printObjectivesToFile("./draw" + "/CMPSODMO" + "/approvePF/" + problem_.getName() + '_' + iteration + ".csv");
			}
			problem_.dynamicChange(iteration);
			dynamicDetected();
			updatePopulation();
			for (int i = 0; i < numberOfSwarm; i++) {
				updatePersonalBest(i);
				updateGlobalBest(i);
			}
			iteration++;
		}
		return Archive;
	} // execute

	private void dynamicDetected() throws JMException, ClassNotFoundException {
		double[][] lastArchive = Archive.writeObjectivesToMatrix();
		for (int i = 0; i < Archive.size(); i++) {
			problem_.evaluate(Archive.get(i));
		}
		double[][] ArchiveObjective = Archive.writeObjectivesToMatrix();
		int times = 0;
		double sums = 0;
		boolean flagChange = false;
		while (times < 0.1 * archiveSize) {
			int random = PseudoRandom.randInt(0, Archive.size() - 1);
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {

				if (Math.abs(lastArchive[random][i] - ArchiveObjective[random][i]) > threshold) {
					System.out.println(iteration);
					SolutionSet temp = new SolutionSet(archiveSize);
					temp = temp.union(Archive);
					Archive.clear();
					for (int j = 0; j < temp.size(); j++) {
						Archive.add(temp.get(j));
					}
					dynamicResponse();
					flagChange = true;
					break;
				}
			}
			if (flagChange) {
				times = archiveSize;
			}
			times++;
		}
	}


	private void dynamicResponse() throws JMException, ClassNotFoundException {
		for (int i = 0; i < numberOfSwarm; i++) {
			int times = 0;
			while (times < 0.2 * particlesSize) {
				int random = PseudoRandom.randInt(0, particlesSize - 1);
				Solution newSolution = new Solution(problem_);
				problem_.evaluate(newSolution);
				particles[i].replace(random, new Solution(newSolution));
				times++;
			}

			for (int j = 0; j < particlesSize; j++) {
				problem_.evaluate(particles[i].get(j));
				personalBest[i].replace(j, new Solution(particles[i].get(j)));
				Archive.add(new Solution(particles[i].get(j)));
			}
			RealMatrix matrix = new Array2DRowRealMatrix(particles[i].writeObjectivesToMatrix());
			int minIndex = matrix.getColumnVector(i).getMinIndex();
			globalBest.replace(i, new Solution(particles[i].get(minIndex)));
		}
	}

	private void updatePersonalBest(int index) {
		for (int i = 0; i < particlesSize; i++) {
			if (personalBest[index].get(i).getObjective(index) >= particles[index].get(i).getObjective(index)) {
				personalBest[index].replace(i, new Solution(particles[index].get(i)));
			}
		}

	}

	private void updateGlobalBest(int index) {
		double min = Double.POSITIVE_INFINITY;
		int minIndex = 0;
		for (int i = 0; i < particles[index].size(); i++) {
			if (particles[index].get(i).getObjective(index) < min) {
				min = particles[index].get(i).getObjective(index);
				minIndex = i;
			}
		}
		globalBest.replace(index, new Solution(particles[index].get(minIndex)));
	}


	private void updatePopulation() throws JMException {
		Variable[] pBest, gBest;
		double r1, r2, r3, r;
		for (int swarmIndex = 0; swarmIndex < numberOfSwarm; swarmIndex++) {
			for (int i = 0; i < particlesSize; i++) {
				Variable[] particle = particles[swarmIndex].get(i).getDecisionVariables();
				double[] velocity = particles[swarmIndex].get(i).getSpeed();
				pBest = personalBest[swarmIndex].get(i).getDecisionVariables();
				int random = PseudoRandom.randInt(0, Archive.size() - 1);
				gBest = globalBest.get(swarmIndex).getDecisionVariables();
				for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
					r1 = PseudoRandom.randDouble();
					r2 = PseudoRandom.randDouble();
					r3 = PseudoRandom.randDouble();
					r = r1 + r2 + r3;
					double w = wInitial - (wInitial - wFinal) * (double) (iteration % 10) / 10;
					double temp = (w * velocity[j]) + c1 * r1 * (pBest[j].getValue() - particle[j].getValue()) / r
							+ c2 * r2 * (gBest[j].getValue() - particle[j].getValue()) / r
							+ c3 * r3 * (Archive.get(random).getDecisionVariables()[j].getValue() - particle[j].getValue()) / r;
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
					if (particles[swarmIndex].get(n).getSpeed(var) > problem_.getUpperLimit(var)) {
						particles[swarmIndex].get(n).setSpeed(var, problem_.getUpperLimit(var));
					}
					if (particles[swarmIndex].get(n).getSpeed(var) < -problem_.getUpperLimit(var)) {
						particles[swarmIndex].get(n).setSpeed(var, -problem_.getUpperLimit(var));
					}

				}
				problem_.evaluate(particles[swarmIndex].get(n));
				Archive.add(new Solution(particles[swarmIndex].get(n)));
			}
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
