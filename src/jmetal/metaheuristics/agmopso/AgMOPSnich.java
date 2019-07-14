package jmetal.metaheuristics.agmopso;
import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;

import java.util.Vector;

public class AgMOPSnich extends Algorithm {
	private static final long serialVersionUID = 2107684627645440737L;
	private Problem problem;
	private int run;
	private double max_d = 0.0;
	private int T_;
	private int[][] neighborhood_;
	private int populationSize;
	//Stores the population
	private SolutionSet leader_ind, cpopulation;
	//Z vector (ideal point)
	double[] idealPoint;
	//Lambda vectors
	double[][] lamdaVectors;
	//Stores the none dominated leaders
	private CrowdingArchive archive;
	Solution[] indArray_;
	// select the aggregation function to be used
	String functionType_;
	Operator crossoverOperator, mutationOperator, cloneoperator;
	// store the number of the particles' evaluations
	int maxIterations;
	private Distance distance_;

	public AgMOPSnich(Problem problem) {
		super(problem);
		this.problem = problem;
	} // MOPSOD

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		functionType_ = "_NBI";//((String) getInputParameter("functionType"));
		maxIterations = ((Integer) this.getInputParameter("maxIterations")).intValue();
		populationSize = ((Integer) this.getInputParameter("swarmSize")).intValue();
		archive = new CrowdingArchive(populationSize, problem.getNumberOfObjectives());
		int clonesize = (int) populationSize / 5;
		T_ = (int) populationSize / 5;
		SolutionSet clonepopulation = new SolutionSet(clonesize);
		int evelations = 0;
		int max_evelations = populationSize * maxIterations;

		cpopulation = new SolutionSet(populationSize);
//		temppopulation=new SolutionSet(populationSize);
		indArray_ = new Solution[problem.getNumberOfObjectives()];
		distance_ = new Distance();
		cloneoperator = operators_.get("clone");
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		neighborhood_ = new int[populationSize][T_];
		idealPoint = new double[problem.getNumberOfObjectives()];
		lamdaVectors = new double[populationSize][problem.getNumberOfObjectives()];
		leader_ind = new SolutionSet(populationSize);

		initUniformWeight();
		initNeighborhood();

		// initialize population
		initPopulation();

		// initialize the Ideal Point
		initIdealPoint(leader_ind);
		find_leader();
		SolutionSet leader_copy = new SolutionSet(populationSize);
		for (int i = 0; i < populationSize; i++) {
			leader_copy.add(new Solution(leader_ind.get(i)));
		}
		leader_copy.sort(new jmetal.util.comparators.CrowdingComparator());
		//get the clone population from the first front
		for (int k = 0; k < leader_copy.size() && k < clonesize; k++) {
			clonepopulation.add(leader_copy.get(k));
		} // for
		// STEP 2. Update
		while (evelations < max_evelations) {
			//1.CLONE POPULATION
			cpopulation = (SolutionSet) cloneoperator.execute(clonepopulation);

			for (int i = 0; i < cpopulation.size(); i++) {
				Solution[] particle2 = new Solution[2];
				int ran;
				particle2[0] = cpopulation.get(i);
				ran = PseudoRandom.randInt(0, cpopulation.size() - 1);
				particle2[1] = cpopulation.get(ran);

				Solution[] offSpring = (Solution[]) crossoverOperator.execute(particle2);
				mutationOperator.execute(offSpring[0]);
				problem.evaluate(offSpring[0]);
				if (problem.getNumberOfConstraints() != 0) {
					problem.evaluateConstraints(offSpring[0]);
				}
				updateReference(offSpring[0]);
				//offSpring[0].setsearch_type(1);
				archive.add(new Solution(offSpring[0]));
				evelations++;
			}
//			archive.printObjectivesToFile(problem_.getName()+"_archive_it"+it);
//			population.printObjectivesToFile(problem_.getName()+"_populaiton_it"+it);
			find_leader();
			//use leaders to guided the population search
			this.ArchiveGuidedSearch();
//			archive.Suppress();
//			distance_.crowdingDistanceAssignment(archive,problem.getNumberOfObjectives());
			find_leader();
			leader_copy = new SolutionSet(populationSize);
			for (int i = 0; i < populationSize; i++) {
				leader_copy.add(new Solution(leader_ind.get(i)));
			}
			leader_copy.sort(new jmetal.util.comparators.CrowdingComparator());
			//get the clone population from the first front
			clonepopulation.clear();
			for (int k = 0; k < leader_copy.size() && k < clonesize; k++) {
				clonepopulation.add(leader_copy.get(k));
			} // for
			evelations += populationSize;
//			archive.printObjectivesToFile(problem.getName()+"_POPit_"+iteration);
		}
		return archive;
	}

	public void initUniformWeight() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			for (int n = 0; n < populationSize; n++) {
				double a = 1.0 * n / (populationSize - 1);
				lamdaVectors[n][0] = a;
				lamdaVectors[n][1] = 1 - a;
				nw++;
			} // for
		} // if
		else if (problem_.getNumberOfObjectives() == 3) {
			int H_ = 13;
			int i, j;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						lamdaVectors[nw][0] = (double) (1.0 * i) / H_;
						lamdaVectors[nw][1] = (double) (1.0 * j) / H_;
						lamdaVectors[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
						nw++;
					} // if
				} // for
			} // for
		} // else
		else if (problem_.getNumberOfObjectives() == 5) {
			int H_ = 6;
			int a, b, c, d;
			for (a = 0; a <= H_; a++) {
				for (b = 0; b <= H_; b++) {
					for (c = 0; c <= H_; c++) {
						for (d = 0; d <= H_; d++) {
							if (a + b + c + d <= H_) {
								lamdaVectors[nw][0] = (double) (1.0 * a) / H_;
								lamdaVectors[nw][1] = (double) (1.0 * b) / H_;
								lamdaVectors[nw][2] = (double) (1.0 * c) / H_;
								lamdaVectors[nw][3] = (double) (1.0 * d) / H_;
								lamdaVectors[nw][4] = (double) (1.0 * (H_ - a - b - c - d) / H_);
								nw++;
							}
						}
					}
				}
			}
		} else if (problem_.getNumberOfObjectives() == 8) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[120][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[36][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										if (a + b + c + d + e + f + g <= H1_) {
											lambda1[nw1][0] = (double) (1.0 * a) / H1_;
											lambda1[nw1][1] = (double) (1.0 * b) / H1_;
											lambda1[nw1][2] = (double) (1.0 * c) / H1_;
											lambda1[nw1][3] = (double) (1.0 * d) / H1_;
											lambda1[nw1][4] = (double) (1.0 * e) / H1_;
											lambda1[nw1][5] = (double) (1.0 * f) / H1_;
											lambda1[nw1][6] = (double) (1.0 * g) / H1_;
											lambda1[nw1][7] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g) / H1_);
											nw1++;
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										if (a + b + c + d + e + f + g <= H2_) {
											lambda2[nw2][0] = (double) (1.0 * a) / H2_;
											lambda2[nw2][1] = (double) (1.0 * b) / H2_;
											lambda2[nw2][2] = (double) (1.0 * c) / H2_;
											lambda2[nw2][3] = (double) (1.0 * d) / H2_;
											lambda2[nw2][4] = (double) (1.0 * e) / H2_;
											lambda2[nw2][5] = (double) (1.0 * f) / H2_;
											lambda2[nw2][6] = (double) (1.0 * g) / H2_;
											lambda2[nw2][7] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g) / H2_);
											nw2++;
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (int i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (int i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		} else if (problem_.getNumberOfObjectives() == 10) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[220][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[55][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g, h, i;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										for (h = 0; h <= H1_; h++) {
											for (i = 0; i <= H1_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H1_) {
													lambda1[nw1][0] = (double) (1.0 * a) / H1_;
													lambda1[nw1][1] = (double) (1.0 * b) / H1_;
													lambda1[nw1][2] = (double) (1.0 * c) / H1_;
													lambda1[nw1][3] = (double) (1.0 * d) / H1_;
													lambda1[nw1][4] = (double) (1.0 * e) / H1_;
													lambda1[nw1][5] = (double) (1.0 * f) / H1_;
													lambda1[nw1][6] = (double) (1.0 * g) / H1_;
													lambda1[nw1][7] = (double) (1.0 * h) / H1_;
													lambda1[nw1][8] = (double) (1.0 * i) / H1_;
													lambda1[nw1][9] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g - h - i) / H1_);
													nw1++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										for (h = 0; h <= H2_; h++) {
											for (i = 0; i <= H2_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H2_) {
													lambda1[nw2][0] = (double) (1.0 * a) / H2_;
													lambda1[nw2][1] = (double) (1.0 * b) / H2_;
													lambda1[nw2][2] = (double) (1.0 * c) / H2_;
													lambda1[nw2][3] = (double) (1.0 * d) / H2_;
													lambda1[nw2][4] = (double) (1.0 * e) / H2_;
													lambda1[nw2][5] = (double) (1.0 * f) / H2_;
													lambda1[nw2][6] = (double) (1.0 * g) / H2_;
													lambda1[nw2][7] = (double) (1.0 * h) / H2_;
													lambda1[nw2][8] = (double) (1.0 * i) / H2_;
													lambda1[nw2][9] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g - h - i) / H2_);
													nw2++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		}
//		for (int i=0;i<nw;i++){
//			for(int j=0;j<problem_.getNumberOfObjectives();j++){
//				if(lambdaVectors[i][j] == 0)
//					lambdaVectors[i][j] = 0.000001;
//			}
//		}
		if (nw != populationSize) {
			System.out.println(nw + "---" + (populationSize));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
		//Apply the WS-transformation on the generated weight vectors
		for (int i = 0; i < populationSize; i++) {
			double prod = 1.0, sum = 0.0;
			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
				prod = prod * lamdaVectors[i][j];
			}
			if (prod != 0.0) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / lamdaVectors[i][j];
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / lamdaVectors[i][j] / sum;
				}
			} else {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / (lamdaVectors[i][j] + 0.0000001);
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / (lamdaVectors[i][j] + 0.0000001) / sum;
				}
			}
		}
	} // initUniformWeight

	public void initNeighborhood() {
		double[] x = new double[populationSize];
		int[] idx = new int[populationSize];

		for (int i = 0; i < populationSize; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize; j++) {
				x[j] = Utils.distVector(lamdaVectors[i], lamdaVectors[j]);
				idx[j] = j;
			} // for
			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize, T_);
			//minfastsort(x,idx,population.size(),niche);
			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize; i++) {
			Solution newSolution = new Solution(problem);
			problem.evaluate(newSolution);
			if (this.problem.getNumberOfConstraints() != 0) {
				problem.evaluateConstraints(newSolution);
			}
			double d = diversity(newSolution, lamdaVectors[i]);
			if (d > max_d) {
				newSolution.setCrowdingDistance(d / max_d);
			}
			leader_ind.add(newSolution);
			archive.add(newSolution);
		}
	} // initPopulation

//	public void orderPopulation(SolutionSet pop) {
//		population = new SolutionSet(populationSize);
//
//		double fitnesses[][] = new double[this.populationSize][this.populationSize];
//		for (int i = 0; i < this.populationSize; i++)
//			for (int j = 0; j < this.populationSize; j++)
//				fitnesses[i][j] = this.fitnessFunction(pop.get(i),
//						this.lambdaVectors[j]);
//		for (int i = 0; i < this.populationSize; i++) {
//			double minFit = Double.MAX_VALUE;
//			int particleIndex = -1;
//			for (int j = 0; j < this.populationSize; j++) {
//				if (fitnesses[j][i] < minFit) {
//					minFit = fitnesses[j][i];
//					particleIndex = j;
//				}
//			}
//			this.population.add(pop.get(particleIndex));
//			for (int n = 0; n < this.populationSize; n++)
//				fitnesses[particleIndex][n] = Double.MAX_VALUE;
//			fitnesses[particleIndex][i] = Double.MAX_VALUE;
//		}
//	}

	public void CD() {
		int best_ind;
		double minFit, fitnesse;
		for (int i = 0; i < archive.size(); i++) {
			best_ind = -1;
			minFit = Double.MAX_VALUE;
			for (int j = 0; j < populationSize; j++) {
				fitnesse = this.fitnessFunction(archive.get(i), this.lamdaVectors[j]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					best_ind = j;//find the best subproblem for archive(i)
				}
			}
			archive.get(i).setCrowdingDistance(diversity(archive.get(i), this.lamdaVectors[best_ind]) / max_d);
		}
	}

	public void find_leader() {
		int best_ind;
		double minFit, fitnesse;
		for (int i = 0; i < this.populationSize; i++) {
			best_ind = -1;
			minFit = Double.MAX_VALUE;
			for (int j = 0; j < archive.size(); j++) {
				fitnesse = this.fitnessFunction(archive.get(j), this.lamdaVectors[i]);
				if (fitnesse < minFit) {
					minFit = fitnesse;
					best_ind = j;
				}
			}
			if (fitnessFunction(leader_ind.get(i), lamdaVectors[i]) > minFit) {
				double d = diversity(archive.get(best_ind), this.lamdaVectors[i]);
				if (d > max_d) {
					;
				}
				max_d = d;
				archive.get(best_ind).setCrowdingDistance(d / max_d);
				leader_ind.replace(i, archive.get(best_ind));
			}
		}
	}

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
				p = PseudoRandom.randInt(0, populationSize - 1);
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

	private void ArchiveGuidedSearch() throws JMException {
		SolutionSet offspringpop = new SolutionSet(populationSize);
		int g1;
		double sign = 1.0;
		Variable[] pbest, lbest1, gbest1;
		for (int n = 0; n < populationSize; n++) {
			//pbest
			pbest = leader_ind.get(n).getDecisionVariables();
			//gbest
			g1 = PseudoRandom.randInt(0, this.archive.size() - 1); // select random  leader
			gbest1 = archive.get(g1).getDecisionVariables();
			//lbest
			Vector<Integer> p = new Vector<Integer>();
			matingSelection(p, n, 1, 1); //Select a neighborhood sub-problem of n
			lbest1 = leader_ind.get(p.get(0)).getDecisionVariables();
			if (fitnessFunction(archive.get(g1), lamdaVectors[n]) < fitnessFunction(leader_ind.get(p.get(0)), lamdaVectors[n])) {
				sign = -1.0;
			}
			double f = diversity(leader_ind.get(n), lamdaVectors[n]) / max_d;
			Solution offspring = new Solution(leader_ind.get(n));
			for (int var = 0; var < problem.getNumberOfVariables(); var++) {
				double temp = pbest[var].getValue() + f * sign * (lbest1[var].getValue() - gbest1[var].getValue());
				if (temp < problem.getLowerLimit(var)) {
					temp = problem.getLowerLimit(var);
				}
				if (temp > problem.getUpperLimit(var)) {
					temp = problem.getUpperLimit(var);
				}
				offspring.getDecisionVariables()[var].setValue(temp);
			}// end for
			offspringpop.add(offspring);
		}
		for (int i = 0; i < this.populationSize; i++) {
			Solution particle = offspringpop.get(i);
			// evaluate the new version of the population
			problem.evaluate(particle);
			if (problem.getNumberOfConstraints() != 0) {
				problem.evaluateConstraints(particle);
			}
			// Update the ideal point
			updateReference(particle);
			// Update of solutions
//			population.replace(i, new Solution(particle)); // change position
			this.archive.add(new Solution(particle));
		}
	}

	void initIdealPoint(SolutionSet pop) throws JMException,
			ClassNotFoundException {
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			idealPoint[i] = 1.0e+30;
			indArray_[i] = new Solution(problem);
			problem.evaluate(indArray_[i]);
			// evaluations++;
		} // for

		for (int i = 0; i < populationSize; i++) {
			updateReference(pop.get(i));
		}

	} // initIdealPoint

	void updateReference(Solution individual) {
		for (int n = 0; n < problem.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < idealPoint[n]) {
				idealPoint[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference

	double fitnessFunction(Solution individual, double[] lamda) {

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;
			for (int n = 0; n < problem.getNumberOfObjectives(); n++) {
				// double diff = Math.abs(individual.getObjective(n)
				// - this.idealPoint[n]);

				double diff = Math.abs(individual.getObjective(n)
						- idealPoint[n]);

				double feval;
				if (lamda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lamda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			return maxFun;
		} // if

		else if (functionType_.equals("_WSUM")) {

			double sum = 0;
			for (int n = 0; n < problem.getNumberOfObjectives(); n++) {
				sum += (lamda[n]) * individual.getObjective(n);
			}
			return sum;

		} // if
		else if ("_NBI".equals(functionType_)) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

			d1 = d2 = nl = 0.0;
			for (i = 0; i < problem.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - idealPoint[i]) * lamda[i];
				nl += Math.pow(lamda[i], 2.0);
			}
			d1 = Math.abs(d1) / Math.sqrt(nl);
			if (nl == 0.0) {
				System.out
						.println("ERROR: dived by zero(bad weihgted vector)\n");
				System.exit(0);
			}
			for (i = 0; i < problem.getNumberOfObjectives(); i++) {
				d2 += Math.pow((individual.getObjective(i) - idealPoint[i])
						- (d1 * lamda[i]), 2.0);
			}
			d2 = Math.sqrt(d2);
			fin = (d1 + theta * d2);
			return fin;
		} else {
			System.out.println("SDMOPSO.fitnessFunction: unknown type "
					+ functionType_);
			return 0;
		}
	} // fitnessEvaluation

	public String getname() {
		return "AGMOPSOnich";
	}
	double diversity(Solution individual, double[] lamda) {
		int i;
		double d1, nl;

		d1 = nl = 0.0;
		for (i = 0; i < problem.getNumberOfObjectives(); i++) {
			d1 += (individual.getObjective(i) - idealPoint[i]) * lamda[i];
			nl += Math.pow(lamda[i], 2.0);
		}
		d1 = Math.abs(d1) / Math.sqrt(nl);
		if (nl == 0.0) {
			System.out
					.println("ERROR: dived by zero(bad weihgted vector)\n");
			System.exit(0);
		}
		return d1;
	} // fitnessEvaluation
} // MOPSOD