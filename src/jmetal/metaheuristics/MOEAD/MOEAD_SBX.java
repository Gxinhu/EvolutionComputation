//  MOEAD.java
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

package jmetal.metaheuristics.MOEAD;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;

import java.util.Vector;

public class MOEAD_SBX extends Algorithm {

	private int populationSize_;
	/**
	 * Stores the population
	 */
	private SolutionSet population_;
	/**
	 * Z vector (ideal point)
	 */
	double[] z_;
	/**
	 * Lambda vectors
	 */
	//Vector<Vector<Double>> lambda_ ;
	double[][] lambda_;
	/**
	 * T: neighbour size
	 */
	int T_;
	int H_ = 12;
	/**
	 * Neighborhood
	 */
	int[][] neighborhood_;
	/**
	 * delta: probability that parent solutions are selected from neighbourhood
	 */
	double delta_;
	/**
	 * nr: maximal number of solutions replaced by each child solution
	 */
	int nr_;
	Solution[] indArray_;
	String functionType_;
	int evaluations_;
	/* if only one layer is adopted, div2_=0 */
	int div1_;  // divisions in the boundary layer
	int div2_;  // divisions in the inside layer
	/**
	 * Operators
	 */
	Operator crossover_;
	Operator mutation_;

	String dataDirectory_;
	int r;

	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public MOEAD_SBX(Problem problem, int run) {
		super(problem);
		r = run;
		functionType_ = "_PBI";

	} // DMOEA

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxEvaluations;

		evaluations_ = 0;
		maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();

		div1_ = ((Integer) this.getInputParameter("div1")).intValue();

		div2_ = ((Integer) this.getInputParameter("div2")).intValue();

		/* generate two-layer weight vectors */
		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();

		/*the population size is the same with the number of weight vectors*/
		populationSize_ = vg.getVectors().length;

//    populationSize = ((Integer) this.getInputParameter("populationSize")).intValue();
//    dataDirectory_ = this.getInputParameter("dataDirectory").toString();
//    System.out.println("POPSIZE: "+ populationSize) ;

		population_ = new SolutionSet(populationSize_);
		indArray_ = new Solution[problem_.getNumberOfObjectives()];

		T_ = 20;
		delta_ = 0.9;
		nr_ = 2;

		neighborhood_ = new int[populationSize_][T_];

		z_ = new double[problem_.getNumberOfObjectives()];
		//lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
//    lambda_ = new double[populationSize][problem_.getNumberOfObjectives()];

		crossover_ = operators_.get("crossover"); // default: DE crossover
		mutation_ = operators_.get("mutation");  // default: polynomial mutation

		// STEP 1. Initialization
		initNeighborhood();

		// STEP 1.2. Initialize population
		initPopulation();

		// STEP 1.3. Initialize z_
		initIdealPoint();
		// STEP 2. Update
		do {
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);

			for (int i = 0; i < populationSize_; i++) {
				int n = permutation[i];
				int type;
				double rnd = PseudoRandom.randDouble();

				// STEP 2.1. Mating selection based on probability
				if (rnd < 0.9) // if (rnd < realb)
				{
					type = 1;   // neighborhood
				} else {
					type = 2;   // whole population
				}
				Vector<Integer> p = new Vector<Integer>();
				matingSelection(p, n, 1, type);//NOTE! We use local search here.

				// STEP 2.2. Reproduction
				Solution child;
				Solution[] parents = new Solution[2];
				parents[1] = population_.get(p.get(0));
				parents[0] = population_.get(n);
				// Apply SBX crossover
				Solution[] offSpring = (Solution[]) crossover_.execute(parents);
				child = offSpring[0];
				// Apply mutation
				mutation_.execute(child);

				// Evaluation
				problem_.evaluate(child);

				evaluations_++;

				// STEP 2.3. Repair. Not necessary

				// STEP 2.4. Update z_
				updateReference(child);

				// STEP 2.5. Update of solutions
				updateProblem(child, n, type);//NOTE: We use global replacement here.
			} // for
		} while (evaluations_ < maxEvaluations);

		return population_;
	}

	/**
	 *
	 */
	public void initNeighborhood() {
		double[] x = new double[populationSize_];
		int[] idx = new int[populationSize_];

		for (int i = 0; i < populationSize_; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize_; j++) {
				x[j] = Utils.distVector(lambda_[i], lambda_[j]);
				//x[j] = dist_vector(population[i].namda,population[j].namda);
				idx[j] = j;
				//System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]: "+idx[j]) ;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize_, T_);
			//minfastsort(x,idx,population.size(),niche);

			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	/**
	 *
	 */
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			evaluations_++;
			population_.add(newSolution);
		} // for
	} // initPopulation

	/**
	 *
	 */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			z_[i] = 1.0e+30;
			indArray_[i] = new Solution(problem_);
			problem_.evaluate(indArray_[i]);
			evaluations_++;
		} // for

		for (int i = 0; i < populationSize_; i++) {
			updateReference(population_.get(i));
		} // for
	} // initIdealPoint

	/**
	 *
	 */
	public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
		// list : the set of the indexes of selected mating parents
		// cid  : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
				//p = population[cid].table[r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			//if (flag) list.push_back(p);
			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection

	/**
	 * @param individual
	 */
	void updateReference(Solution individual) {
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < z_[n]) {
				z_[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference

	/**
	 * @param id
	 * @param type
	 */
	void updateProblem(Solution indiv, int id, int type) {
		// indiv: child solution
		// id:   the id of current subproblem
		// type: update solutions in - neighborhood (1) or whole population (otherwise)
		int size;
		int time;

		time = 0;

		if (type == 1) {
			size = neighborhood_[id].length;
		} else {
			size = population_.size();
		}
		int[] perm = new int[size];

		Utils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[id][perm[i]];
			} else {
				k = perm[i];      // calculate the values of objective function regarding the current subproblem
			}
			double f1, f2;

			f1 = fitnessFunction(population_.get(k), lambda_[k]);
			f2 = fitnessFunction(indiv, lambda_[k]);

			if (f2 < f1) {
				population_.replace(k, new Solution(indiv));
				//population[k].indiv = indiv;
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed 'limit'
			if (time >= nr_) {
				return;
			}
		}
	} // updateProblem

	/**
	 * Calculate the dot product of two vectors
	 *
	 * @param vec1
	 * @param vec2
	 * @return
	 */
	public double innerproduct(double[] vec1, double[] vec2) {
		double sum = 0;

		for (int i = 0; i < vec1.length; i++) {
			sum += vec1[i] * vec2[i];
		}

		return sum;
	}

	/**
	 * Calculate the norm of the vector
	 *
	 * @param z
	 * @return
	 */
	public double norm_vector(double[] z) {
		double sum = 0;

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			sum += z[i] * z[i];
		}

		return Math.sqrt(sum);
	}

	double fitnessFunction(Solution indiv, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(indiv.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} else if (functionType_.equals("_TCHE2")) {
			double maxFun = -1.0e+30;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				double diff = Math.abs(indiv.getObjective(i) - z_[i]);

				double feval;
				if (lambda[i] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[i];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			fitness = maxFun;
		} else if (functionType_.equals("_PBI")) {
			double theta; // penalty parameter
			theta = 5.0;

			// normalize the weight vector (line segment)
			double nd = norm_vector(lambda);
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				lambda[i] = lambda[i] / nd;
			}

			double[] realA = new double[problem_.getNumberOfObjectives()];
			double[] realB = new double[problem_.getNumberOfObjectives()];

			// difference between current point and reference point
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				realA[n] = (indiv.getObjective(n) - z_[n]);
			}

			// distance along the line segment
			double d1 = Math.abs(innerproduct(realA, lambda));

			// distance to the line segment
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				realB[n] = (indiv.getObjective(n) - (z_[n] + d1 * lambda[n]));
			}
			double d2 = norm_vector(realB);

			fitness = d1 + theta * d2;
		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation
} // MOEAD

