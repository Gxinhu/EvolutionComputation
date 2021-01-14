/**
 * PaRPEA.java
 * In this class, the admoea is combined with MAF as in VaEA
 * This is an MOEA based on adoptive shape method, a release version
 * Created 6-4-2017
 * Last modified 8-7-2017
 * Code by Dr. Xiang (gzhuxiang_yi@163.com)
 */
package jmetal.metaheuristics.parpea;

import java.util.Comparator;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.myutils.CaptureConvergence;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Utils;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.FitnessComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import Jama.Matrix;

public class PaRPEA extends Algorithm {

    private int populationSize_;
    private int objectives_;
    private double q_; // Curvature
    private double eps_ = 0.05; //

    private SolutionSet population_;
    private SolutionSet offspringPopulation_;
    private SolutionSet union_;

    private double[] zp_; // Normalized ideal point
    private double[] nzp_; // Normalized nadir point
    private double[] referencePoint_; // Reference point

    //-----------Begin:The following four variables are used in normalization---------
    double[] zideal; // Original ideal point
    double[] zmax;   // Original nadir point
    double[][] extremePoints;
    double[] intercepts;
    //-----------------------------------End------------------------------------

    /*------- Begin:The following three variables are used in the environmental selection	 */
    double[][] AngleMatrix;
    boolean[] removed;
    double[] minAngleArray;
    /*-----------------------------End:-------------------*/

    int[] rpTypeCounter = new int[3]; // rpTypeCounter[0]--convex; rpTypeCounter[1]--linear; rpTypeCounter[2]--concave
    // Operators
    Operator crossover_;
    Operator mutation_;
    Operator selection_;

    int evaluations_;
    int maxEvaluations_;
    int iteration_;
    int maxIteration_;

    boolean normalize_; // Do normalization or not
    CaptureConvergence captureConvergence_; // our class
    boolean convergenceFlag = false;

    /**
     * Constructor
     *
     * @param problem
     */
    public PaRPEA(Problem problem) {
        super(problem);
        objectives_ = problem.getNumberOfObjectives();
        zideal = new double[objectives_];
        zmax = new double[objectives_];
        zp_ = new double[objectives_];
        nzp_ = new double[objectives_];
        referencePoint_ = new double[objectives_];

        if (convergenceFlag == true)
            captureConvergence_ = new CaptureConvergence(problem.getName(), "PaRPEA+eliminate",
                    problem.getNumberOfObjectives());
    } // PaRPEA

    public SolutionSet execute() throws JMException, ClassNotFoundException {
        evaluations_ = 0;
        iteration_ = 0;
        populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
        int maxGeneration = ((Integer) this.getInputParameter("maxGenerations")).intValue();
        maxEvaluations_ = maxGeneration * populationSize_;
        normalize_ = ((Boolean) this.getInputParameter("normalize")).booleanValue();

        mutation_ = operators_.get("mutation");
        crossover_ = operators_.get("crossover");
        selection_ = operators_.get("selection");

        maxIteration_ = maxEvaluations_ / populationSize_;

        /**
         * Step 1: Initialization of a random population
         */
        initPopulationRandom();

        if (convergenceFlag == true)
            captureConvergence_.runCaptureConvergence(0, population_);

        // Main loop
        while (evaluations_ <= maxEvaluations_) {

            /**
             * Step 2. Create an offspring population
             */
            createOffspringPopulation();
            /**
             * Step 3: Environmental selection
             */
            environmentalSelection();

            if (convergenceFlag == true) {
                captureConvergence_.runCaptureConvergence(evaluations_, population_);
            }

            iteration_++;

            if (evaluations_ >= maxEvaluations_)
                break;

        } // while

        System.out.println("Convex = " + rpTypeCounter[0] + "  ,Linear = " + rpTypeCounter[1] + ",Concave = " + rpTypeCounter[2]);
        System.out.println("PaRPEA popSize = " + population_.size());
        System.out.println("evaluations = " + evaluations_);
        Ranking ranking = new NondominatedRanking(population_);

        // Return non-dominated solutions
        SolutionSet nondominatedSet = ranking.getSubfront(0);
        System.out.println("The number of nondominated solutions = " + nondominatedSet.size());

        return nondominatedSet;

    } // execute

    /**
     * Generate initial population randomly
     *
     * @throws JMException
     * @throws ClassNotFoundException
     */
    public void initPopulationRandom() throws JMException, ClassNotFoundException {

        /**
         * Method 1: Random generate N solutions
         */
        population_ = new SolutionSet(populationSize_);

        for (int i = 0; i < populationSize_; i++) {
            Solution newSolution = new Solution(problem_);

            problem_.evaluate(newSolution);
            problem_.evaluateConstraints(newSolution);
//			evaluations_++;

            population_.add(newSolution);
        } // for

    } // initPopulationRandom

    /**
     * For each solution i, generate 2 new solutions, and stored them in offspringPopulation_
     *
     * @throws JMException
     */
    public void createOffspringPopulation() throws JMException {
        offspringPopulation_ = new SolutionSet(populationSize_);
        /**
         *  Method 1: Crossover and mutation
         */
        for (int i = 0; i < populationSize_ / 2; i++) {

            Solution[] parents = new Solution[2];
            parents[0] = population_.get(i);

            int secParent = 0;
            do {
                secParent = PseudoRandom.randInt(0, populationSize_ - 1);
            } while (secParent == i);

            parents[1] = population_.get(secParent);

//			 parents[1] = (Solution) selection_.execute(new Object[] { i , population_});			
//			 parents[0] = (Solution) selection_.execute(population_);
//		     parents[1] = (Solution) selection_.execute(population_);

            Solution[] offSpring = (Solution[]) crossover_.execute(parents);

            mutation_.execute(offSpring[0]);
            mutation_.execute(offSpring[1]);

            problem_.evaluate(offSpring[0]);
            problem_.evaluateConstraints(offSpring[0]);
            evaluations_++;

            problem_.evaluate(offSpring[1]);
            problem_.evaluateConstraints(offSpring[1]);
            evaluations_++;

            offspringPopulation_.add(offSpring[0]);
            offspringPopulation_.add(offSpring[1]);

        } // for i

        /**
         * Method 2: DE Crossover & mutation
         */
//		offspringPopulation_ = new SolutionSet(populationSize_);
//		 for (int i = 0; i < populationSize_; i++) {
//		
//			 Solution[] parents = new Solution[3];
//			 Solution child;
//			
//			 int p1 = PseudoRandom.randInt(0, populationSize_ - 1);
//			 parents[0] = population_.get(p1);
//			
//			 int p2 = PseudoRandom.randInt(0, populationSize_ - 1);
//			
//			 while (p2 == p1) {
//				 p2 = PseudoRandom.randInt(0, populationSize_ - 1);
//			 }
//			
//			 parents[1] = population_.get(p2);
//			
//			 int p3 = PseudoRandom.randInt(0, populationSize_ - 1);
//			
//			 while (p3 == p1 || p3 == p2) {
//				 p3 = PseudoRandom.randInt(0, populationSize_ - 1);
//			 }

//			 parents[2] = population_.get(p3);

//			 parents = (Solution []) selection_.execute(new Object[] {population_, i});	
//			 
//			 // Apply DE crossover
//			 child = (Solution) crossover_.execute(new Object[] { population_.get(i), parents});
//			 // Apply mutation
//			 mutation_.execute(child); // 
//			
//			 problem_.evaluate(child);
//			 problem_.evaluateConstraints(child);
//			 evaluations_++;
//			 offspringPopulation_.add(child);
//		 } // for i
    } // createOffspringPopulation

    /**
     * Do environmental selection procedure
     */
    public void environmentalSelection() {

        // Step 1: Get the merged population
        union_ = population_.union(offspringPopulation_);

        // Step 2: Split out dominated and nondominated solutions
        SolutionSet[] splitResults = splitDominatedNondominated(union_);
        SolutionSet nondominated = splitResults[0];
        SolutionSet dominated = splitResults[1];

        // Step 3: Calculate curvature based on nondominated solutions
        computeCurvature(nondominated); // Compute curvature q_

        // Step 4: Update reference point
        updateReferencePointBasedOnCurvature();

        population_.clear();

        if (nondominated.size() == populationSize_) { // Occur hardly
            // Add all non-dominated solutions into population_
            for (int i = 0; i < nondominated.size(); i++) {
                population_.add(nondominated.get(i));
            }
        } else if (nondominated.size() < populationSize_) {
            // First add all non-dominated solutions into population_
            for (int i = 0; i < nondominated.size(); i++) {
                population_.add(nondominated.get(i));
            }

            /** --Renormalized may be unnecessary, however it will not affect the runtime significantly--*/
            computeIdealPoint(union_);
            computeMaxPoint(union_);
            computeExtremePoints(union_);
            computeIntercepts();
            normalizePopulation(union_);

            calculateFitnessBasedQ(dominated);
            initializeAngleMatrix(dominated);
            eliminateWithoutExtremeSolution(dominated, population_);

        } else {
			
				/* Used in PaRPEAMAFNoRank�� ����ʹ��classficationbydominance
				nondominated = new SolutionSet(union_.size());
				nondominated = union_;
				normalizePopulation(nondominated);
				
				// Compute norm of each solution
				for (int i = 0; i < nondominated.size(); i++) {
					computeDistance2IdealPoint(nondominated.get(i));
				}
				*/

            SolutionSet[] betterAndWorse = splitBetterWorse(nondominated);
            SolutionSet better = betterAndWorse[0];
            SolutionSet worse = betterAndWorse[1];

            if (better.size() <= populationSize_) {

                for (int i = 0; i < better.size(); i++) {
                    population_.add(better.get(i));
                }

                if (population_.size() == populationSize_) return;

                calculateFitnessBasedQ(worse);
                worse.sort(new FitnessComparator());
                int t = 0;
                while (population_.size() < populationSize_) {
                    population_.add(worse.get(t));
                    t++;
                }

            } else {

                calculateFitnessBasedQ(better);
                if (worse.size() > 0) {
                    initializeAngleMatrix(better);
                    eliminate(better, population_);
                } else {
                    MAF(better);
                }

            } // if

        }

    } // environmentalSelection

    public void updateReferencePoint() {
        for (int i = 0; i < this.objectives_; i++) {
            referencePoint_[i] = 1.05;

        }

    }

    /**
     * Split out nondominated and dominated solutions from the union set
     *
     * @param union
     * @return
     */
    public SolutionSet[] splitDominatedNondominated(SolutionSet union) {
        SolutionSet[] results = new SolutionSet[2];

        boolean[] dominated = flagsDominateNondominate(union);
        int noOfNondominated = 0;

        for (int i = 0; i < union.size(); i++) {
            if (dominated[i] == false) {
                noOfNondominated++;
            }
        }

        results[0] = new SolutionSet(noOfNondominated); // results[0] is the set
        // of nondominated

        for (int i = 0; i < union.size(); i++) {
            if (dominated[i] == false) {
                results[0].add(union.get(i));
            }
        }

        int noOfdominated = 0;

        for (int i = 0; i < union.size(); i++) {
            if (dominated[i] == true) {
                noOfdominated++;
            }
        }

        results[1] = new SolutionSet(noOfdominated); // results[1] is the set of
        // dominated

        for (int i = 0; i < union.size(); i++) {
            if (dominated[i] == true) {
                results[1].add(union.get(i));
            }
        }

        return results;
    }

    /**
     * Split out better and worse solutions from the union set
     *
     * @param union
     * @return
     */
    public SolutionSet[] splitBetterWorse(SolutionSet union) {
        SolutionSet[] results = new SolutionSet[2];

        boolean[] flagWorse = new boolean[union.size()];
        int noOfWorse = 0;


        for (int i = 0; i < union.size(); i++) {
            Solution solution = union.get(i);

            for (int j = 0; j < objectives_; j++) {
                double val = 1.0 + eps_;
                if (solution.getTranslatedObjectives(j) > val) {
                    flagWorse[i] = true;
                    noOfWorse++;
                    break;
                }
            }
        }

        results[0] = new SolutionSet(union.size() - noOfWorse); // results[0] is better ones
        results[1] = new SolutionSet(noOfWorse); // results[1] is worse ones

        for (int i = 0; i < union.size(); i++) {
            if (flagWorse[i] == false) {
                results[0].add(union.get(i));
            } else {
                results[1].add(union.get(i));
            }
        }

        return results;
    }

    /**
     * Compute Curvature based on nondominated solutions
     */
    public void computeCurvature(SolutionSet nondominated) {

        // Use solutions in nondominated set to find zmin and zmax
        computeIdealPoint(nondominated); // zideal
        computeMaxPoint(nondominated);  // zmax

        computeExtremePoints(nondominated); // extreme solutions
        computeIntercepts(); //  intercepts
        normalizePopulation(nondominated);

        double[] V = new double[objectives_];
        for (int i = 0; i < objectives_; i++) {
//			 V[i] = 1.0/objectives_;
            V[i] = 1.0;
        }

        // Compute norm of each solution
        for (int i = 0; i < nondominated.size(); i++) {
            computeDistance2IdealPoint(nondominated.get(i));
        }

//		/**
//		 * Method 1, exact one 
//		 */
//		// Find the closest solution to V based on angle
//		double minTheta = 1.0e+30;
//		int minThetaID = -1;
//		double minThetaNorm = 0.0;
//
//		for (int i = 0; i < nondominated.size(); i++) {
//			double angle = calAngleFromIdealPoint(nondominated.get(i), V);			
//			if (angle < minTheta) {
//				minTheta = angle;
//				minThetaID = i;
//				minThetaNorm = nondominated.get(i).getDistanceToIdealPoint();
//			}
//		}
//			
//		q_ = Math.log(1.0 / objectives_)/ Math.log(minThetaNorm / Math.sqrt(objectives_));
        /*--------------------------------------------------------------------------------*/

        /**
         * Method 2, approximated one
         */
        double[] angles = new double[nondominated.size()];
        int idx[] = new int[nondominated.size()];

        for (int i = 0; i < nondominated.size(); i++) {
            angles[i] = calAngleFromIdealPoint(nondominated.get(i), V);
        }

        Utils.minFastSort(angles, idx, nondominated.size(), this.objectives_);
        double meanNorm = 0.0;

        int no = this.objectives_;
        if (nondominated.size() < this.objectives_) {
            no = nondominated.size();
        }

        for (int i = 0; i < no; i++) {
            meanNorm += nondominated.get(idx[i]).getDistanceToIdealPoint();
        }

        meanNorm = meanNorm / no;

        q_ = meanNorm * Math.sqrt(objectives_);
//        q_ = 1.0; // for ConvexDTLZ3

        /*
         *  ��¼������Ʋο��������
         *   rpTypeCounter[0]--convex; rpTypeCounter[1]--linear; rpTypeCounter[2]--concave
         */
        if (q_ < 0.9)
            rpTypeCounter[0]++;
        else if (q_ > 1.1)
            rpTypeCounter[2]++;
        else
            rpTypeCounter[1]++;

    } // computeCurvature

    /**
     * Update reference point
     */
    public void updateReferencePointBasedOnCurvature() {

//		if (q_ >= 0.9 && q_ <= 1.1 )  {
//			referencePoint_ = zp_;			
//		} else if (q_ > 1.1) {
//			referencePoint_ = zp_;				
//		} else {	
//			referencePoint_ = intercepts; // nor nzp_				
//		}	
        for (int i = 0; i < this.objectives_; i++) {
            if (q_ >= 0.9) {
                referencePoint_[i] = -eps_;
            } else {
                referencePoint_[i] = 1.0 + eps_;
            }

        } // for

    } // updateReferencePointBasedOnCurvature

    void computeExtremePoints(SolutionSet solutionSet) {
        extremePoints = new double[objectives_][objectives_];

        for (int j = 0; j < objectives_; j++) {
            int index = -1;
            double min = Double.MAX_VALUE;

            for (int i = 0; i < solutionSet.size(); i++) {
                double asfValue = asfFunction(solutionSet.get(i), j);
                if (asfValue < min) {
                    min = asfValue;
                    index = i;
                }
            }

            for (int k = 0; k < objectives_; k++)
                extremePoints[j][k] = solutionSet.get(index).getObjective(k);
        }
    }

    double asfFunction(Solution sol, int j) {
        double max = Double.MIN_VALUE;
        double epsilon = 1.0E-6;

        for (int i = 0; i < objectives_; i++) {

            double val = Math.abs(sol.getObjective(i) - zideal[i]);

            if (j != i)
                val = val / epsilon;

            if (val > max)
                max = val;
        }

        return max;
    }

    void computeIntercepts() {

        intercepts = new double[objectives_];

        double[][] temp = new double[objectives_][objectives_];

        for (int i = 0; i < objectives_; i++) {
            for (int j = 0; j < objectives_; j++) {
                double val = extremePoints[i][j] - zideal[j];
                temp[i][j] = val;
            }
        }

        Matrix EX = new Matrix(temp);

        if (EX.rank() == EX.getRowDimension()) {
            double[] u = new double[objectives_];
            for (int j = 0; j < objectives_; j++)
                u[j] = 1;

            Matrix UM = new Matrix(u, objectives_);

            Matrix AL = EX.inverse().times(UM);

            int j = 0;
            for (j = 0; j < objectives_; j++) {

                double aj = 1.0 / AL.get(j, 0) + zideal[j];

                if ((aj > zideal[j]) && (!Double.isInfinite(aj))
                        && (!Double.isNaN(aj)))
                    intercepts[j] = aj;
                else
                    break;
            }
            if (j != objectives_) {
                for (int k = 0; k < objectives_; k++)
                    intercepts[k] = zmax[k];
            }

        } else {
            for (int k = 0; k < objectives_; k++)
                intercepts[k] = zmax[k];
        }

    }

    void normalizePopulation(SolutionSet solutionSet) {
        for (int i = 0; i < solutionSet.size(); i++) {
            Solution sol = solutionSet.get(i);

            for (int j = 0; j < objectives_; j++) {

                double val = (sol.getObjective(j) - zideal[j])
                        / (intercepts[j] - zideal[j]);

                sol.setTranslatedObjectives(val, j);
            }// for
        }

        for (int j = 0; j < objectives_; j++) {
            intercepts[j] = 1.0;
        }// for

    }

    public void calculateFitness(SolutionSet union) {

        // System.out.println("linear");
        for (int i = 0; i < union.size(); i++) {
            double sum = 0.0;
            double dis = 0.0;

            Solution sol = union.get(i);

            for (int j = 0; j < objectives_; j++) {
                sum = sum + sol.getTranslatedObjectives(j);
                dis = dis + (sol.getTranslatedObjectives(j) - referencePoint_[j])
                        * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
            }

            dis = Math.sqrt(dis);
            sol.setDistanceToReferencePoint(dis);
            sol.setFitness(sum);
        }
    }

    public void calculateFitnessBasedQ(SolutionSet union) {
        if (q_ >= 0.9 && q_ <= 1.1) { // linear, the sum of all the objectives
            // (sum for short)
            // System.out.println("linear");
            for (int i = 0; i < union.size(); i++) {
                double sum = 0.0;
                double dis = 0.0;

                Solution sol = union.get(i);

                for (int j = 0; j < objectives_; j++) {
                    sum = sum + sol.getTranslatedObjectives(j);
                    dis = dis + (sol.getTranslatedObjectives(j) - referencePoint_[j])
                            * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
                }

                dis = Math.sqrt(dis);
                sol.setDistanceToReferencePoint(dis);
                sol.setFitness(sum);
            }

        } else if (q_ > 1.1) { // concave ��, Euclidean distance to an ideal
            // point (EdI for short)
            // System.out.println("concave��");
            for (int i = 0; i < union.size(); i++) {
                double sum = 0.0;
                Solution sol = union.get(i);

                for (int j = 0; j < objectives_; j++) {
                    sum = sum + (sol.getTranslatedObjectives(j) - referencePoint_[j])
                            * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
                }
                sum = Math.sqrt(sum);
                sol.setDistanceToReferencePoint(sum);
                sol.setFitness(sum);
            }

        } else { // convex͹,Euclidean distance to the Nadir point (EdN for
            // short)
            // System.out.println("convex͹");
//			double maxDis = 0.0;
//			
//			for (int j = 0; j < objectives_; j++) {
//				maxDis = maxDis + nzp_[j] * nzp_[j];					
//			}
//			maxDis = Math.sqrt(maxDis);	

            for (int i = 0; i < union.size(); i++) {
                double sum = 0.0;
                Solution sol = union.get(i);

                for (int j = 0; j < objectives_; j++) {
                    sum = sum + (sol.getTranslatedObjectives(j) - referencePoint_[j])
                            * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
                }

                sum = Math.sqrt(sum);
                sol.setDistanceToReferencePoint(sum);

                sum = 1.0 / sum;
                sol.setFitness(sum);

            }
        }

    }

    /**
     * Initialize angle matrix
     */
    public void initializeAngleMatrix(SolutionSet solutionSet) {
        // solutionSet.sort(new FitnessComparator());
        AngleMatrix = new double[solutionSet.size()][solutionSet.size()];
        for (int i = 0; i < solutionSet.size(); i++) {
            solutionSet.get(i).setID(i);// Give a solution a unique ID
        }

        for (int i = 0; i < solutionSet.size(); i++) {

            Solution soli = solutionSet.get(i);

            for (int j = i; j < solutionSet.size(); j++) {

                if (i == j) {
                    AngleMatrix[i][i] = 0.0;
                } else {
                    Solution solj = solutionSet.get(j);
                    AngleMatrix[i][j] = calAngle(soli, solj);
                    AngleMatrix[j][i] = AngleMatrix[i][j];
                }

            } // for j

        } // for i

    } // initializeAngleMatrix

    /**
     * Initialize angle matrix from the nzp point
     */
    public void initializeAngleMatrixFromNZP(SolutionSet solutionSet) {
        // solutionSet.sort(new FitnessComparator());
        AngleMatrix = new double[solutionSet.size()][solutionSet.size()];
        for (int i = 0; i < solutionSet.size(); i++) {
            solutionSet.get(i).setID(i);// Give a solution a unique ID
        }

        for (int i = 0; i < solutionSet.size(); i++) {

            Solution soli = solutionSet.get(i);

            for (int j = i; j < solutionSet.size(); j++) {

                if (i == j) {
                    AngleMatrix[i][i] = 0.0;
                } else {
                    Solution solj = solutionSet.get(j);
                    AngleMatrix[i][j] = calAngleFromNZP(soli, solj);
                    AngleMatrix[j][i] = AngleMatrix[i][j];
                }

            } // for j

        } // for i

    } // initializeAngleMatrix

    public boolean[] flagsDominateNondominate(SolutionSet solutionSet) {
        // Find dominated members in population_
        boolean[] dominated = new boolean[solutionSet.size()];
        Comparator dominance = new DominanceComparator();
        // Remove dominated solutions in the population
        for (int i = 0; i < solutionSet.size(); i++) {

            if (dominated[i] == false) { //
                Solution sol_i = solutionSet.get(i);

                for (int j = i + 1; j < solutionSet.size(); j++) {
                    if (dominated[j] == true) {
                        continue;
                    }

                    Solution sol_j = solutionSet.get(j);
                    int flagDominate = dominance.compare(sol_i, sol_j);
                    if (flagDominate == -1) { // sol_i dominates sol_j
                        dominated[j] = true;
                    } else if (flagDominate == 1) { // sol_i is dominated by
                        // sol_j
                        dominated[i] = true;
                        break;
                    }
                }
            }
        }
        return dominated;
    }

    /**
     * Eliminate without considering extreme solutions
     *
     * @param union
     * @param pop
     */
    public void eliminateWithoutExtremeSolution(SolutionSet union, SolutionSet pop) {
        union.sort(new FitnessComparator());
        // Step: 1 Initialize angle matrix
        removed = new boolean[union.size()];
        minAngleArray = new double[union.size()];

        boolean removedFlag = false;

        // Find the minimum angle for each unadded solution
        for (int i = 0; i < union.size(); i++) {
            if (removed[i] == false) {
                double minAng = 1.0e+30;
                int minAngID = -1;

                for (int j = 0; j < union.size(); j++) {

                    if (j == i)
                        continue;

                    if (removed[j] == false && AngleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
                        minAng = AngleMatrix[union.get(i).getID()][union.get(j).getID()];
                        minAngID = j;
                    }

                } // for j

                union.get(i).setClusterID(minAngID);
                union.get(i).setAssociateDist(minAng);
                minAngleArray[i] = minAng;

            } // if
        } // for i

        int remain = populationSize_ - pop.size();

        // Remove solutions
        int numberOfLoops = (union.size()) - remain;
        if (removedFlag == true)
            numberOfLoops = (union.size() - this.objectives_) - remain;

        // if (extremeAdded == true)
        // numberOfLoops = (union.size()-this.objectives_) - remain;

        for (int r = 0; r < numberOfLoops; r++) {
            double minValue = 1.0e+30;
            int minValueID = -1;

            for (int i = 0; i < union.size(); i++) {
                if (removed[i] == true)
                    continue;
                if (minAngleArray[i] < minValue) {
                    minValue = minAngleArray[i];
                    minValueID = i;
                }
            } // for i

            int associatedId = union.get(minValueID).getClusterID();
            int removedID = -1;

            if (union.get(minValueID).getFitness() < union.get(associatedId).getFitness()) {

                removedID = associatedId;

            } else if (union.get(minValueID).getFitness() > union.get(associatedId).getFitness()) {

                removedID = minValueID;

            } else {

                if (PseudoRandom.randDouble() < 0.5) {
                    removedID = associatedId;
                } else {
                    removedID = minValueID;
                }


            } // if

            removed[removedID] = true;

            // Update angles
            for (int i = 0; i < union.size(); i++) {
                if (removed[i] == true)
                    continue;

                if (AngleMatrix[union.get(i).getID()][union.get(removedID).getID()] == union.get(i).getAssociateDist()) { // �������׼ȷ��

                    double minAng = 1.0e+30;
                    int minAngID = -1;

                    for (int j = 0; j < union.size(); j++) {

                        if (j == i)
                            continue;

                        if (removed[j] == false
                                && AngleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
                            minAng = AngleMatrix[union.get(i).getID()][union.get(j).getID()];
                            minAngID = j;
                        }

                    } // for j

                    union.get(i).setClusterID(minAngID);
                    union.get(i).setAssociateDist(minAng);
                    minAngleArray[i] = minAng;
                }
            }

        } // for r

        // Add remain solutions into pop
        for (int i = 0; i < union.size(); i++) {
            if (removed[i] == false) { // && pop.size() < populationSize_
                pop.add(union.get(i));
            }
        } // for i

    } // eliminateWithoutExtremeSolution

    /**
     * Delete solutions
     *
     * @param union
     * @param pop
     */
    public void eliminate(SolutionSet union, SolutionSet pop) {
//		System.out.println("eliminate");
        union.sort(new FitnessComparator());
        // Step: 1 Initialize angle matrix
        removed = new boolean[union.size()];
        minAngleArray = new double[union.size()];
        boolean[] isExtreme = new boolean[union.size()];
        boolean[] considered = new boolean[union.size()];
        int[] removedSolutions = new int[populationSize_];
        int noOfRemoved = 0;

        // Add extreme solutions
//	    int [] perm = (new jmetal.util.PermutationUtility()).intPermutation(objectives_);		
        for (int k = 0; k < objectives_; k++) {
            double[] vector = new double[objectives_];
            vector[k] = 1.0;
            double minAngle = 1e+30;
            int minAngleID = -1;

            for (int i = 0; i < union.size(); i++) {

                if (removed[i] == false) {
                    Solution soli = union.get(i);
//					double angle = calAngle(soli, vector);
                    double angle = this.calAngleFromIdealPoint(soli, vector);
                    if (angle < minAngle) {
                        minAngle = angle;
                        minAngleID = i;
                    }
                }
            } // for i
            pop.add(new Solution(union.get(minAngleID)));
            removed[minAngleID] = true;
            isExtreme[minAngleID] = true;
            considered[minAngleID] = true;
        } // for k


        boolean removedFlag = false;

        // Find the minimum angle for each unadded solution
        for (int i = 0; i < union.size(); i++) {
            if (removed[i] == false) {
                double minAng = 1.0e+30;
                int minAngID = -1;

                for (int j = 0; j < union.size(); j++) {

                    if (j == i)
                        continue;

                    if (removed[j] == false
                            && AngleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
                        minAng = AngleMatrix[union.get(i).getID()][union.get(j).getID()];
                        minAngID = j;
                    }

                } // for j

                union.get(i).setClusterID(minAngID);
                union.get(i).setAssociateDist(minAng);
                minAngleArray[i] = minAng;
            } // if
        } // for i

        int remain = populationSize_ - pop.size();
//		System.out.println("remain" + remain);
        // Remove solutions
        int numberOfLoops = union.size() - pop.size() - remain;
//		System.out.println("numberOfLoops" + numberOfLoops);

        if (removedFlag == true)
            numberOfLoops = (union.size() - this.objectives_) - pop.size() - remain;

        for (int r = 0; r < numberOfLoops; r++) {

            double minValue = 1.0e+30;
            int minValueID = -1;

            for (int i = 0; i < union.size(); i++) {
                if (removed[i] == true)
                    continue;
                if (minAngleArray[i] <= minValue) { // ȡС�ں�
                    minValue = minAngleArray[i];
                    minValueID = i;
                }
            } // for i

            int associatedId = union.get(minValueID).getClusterID();
            int removedID = -1;

            if (union.get(minValueID).getFitness() < union.get(associatedId).getFitness()) {

                removedID = associatedId;

            } else if (union.get(minValueID).getFitness() > union.get(associatedId).getFitness()) {

                removedID = minValueID;

            } else {

//				if (union.get(minValueID).getDistanceToIdealPoint() < union.get(associatedId).getDistanceToIdealPoint()){
//					removedID = associatedId;
//				} else if (union.get(minValueID).getDistanceToIdealPoint() > union.get(associatedId).getDistanceToIdealPoint()){
//					removedID = minValueID;					
//				} else {

                if (PseudoRandom.randDouble() < 0.5) {
                    removedID = associatedId;
                } else {
                    removedID = minValueID;
                }

//				}// IF 

            } // if

            removed[removedID] = true;

            removedSolutions[noOfRemoved] = removedID;
            noOfRemoved++;

            considered[minValueID] = true;
            considered[associatedId] = true;

            // Update angles
            for (int i = 0; i < union.size(); i++) {
                if (removed[i] == true)
                    continue;

                if (AngleMatrix[union.get(i).getID()][union.get(removedID).getID()] == union.get(i).getAssociateDist()) { // �������׼ȷ��
                    double minAng = 1.0e+30;
                    int minAngID = -1;

                    for (int j = 0; j < union.size(); j++) {

                        if (j == i)
                            continue;

                        if (removed[j] == false && AngleMatrix[union.get(i).getID()][union.get(j).getID()] < minAng) {
                            minAng = AngleMatrix[union.get(i).getID()][union.get(j).getID()];
                            minAngID = j;
                        }

                    } // for j

//					if (minAngID == -1) {
//						System.out.println("union.size = " + union.size());
//						System.out.println("r = " + r);
//						System.out.println("numberOfLoops = " + numberOfLoops);
//						System.out.println(minAng);
//						System.out.println("here minAngID == -1");
//					}

                    union.get(i).setClusterID(minAngID);
                    union.get(i).setAssociateDist(minAng);
                    minAngleArray[i] = minAng;
                }
            }

        } // for r

        // Add remain solutions into pop
        for (int i = 0; i < union.size(); i++) {
            if (removed[i] == false && isExtreme[i] == false) { //
                pop.add(union.get(i));
            }
        } // for i

    } // niching

    /**
     * Maximum angle first principle
     *
     * @param lastFront
     */
    public void MAF(SolutionSet lastFront) {

//		System.out.println("MAF");
        int n = lastFront.size();
        double[] angles = new double[n];
        int[] index = new int[n];
        boolean[] removed = new boolean[n];
        lastFront.sort(new FitnessComparator());
        int remain = populationSize_;


        for (int o = 0; o < this.objectives_; o++) {
            double minAngle2Axis = 1.0e+30;
            int minAngle2AxisID = -1;

            for (int i = 0; i < n; i++) {
                if (removed[i] == false) {
                    Solution solLastFront = lastFront.get(i);
                    double angle = Math.acos(Math.abs(solLastFront.getTranslatedObjectives(o) / solLastFront.getDistanceToIdealPoint()));

                    if (angle < minAngle2Axis) {
                        minAngle2Axis = angle;
                        minAngle2AxisID = i;
                    }
                }
            }// for

            removed[minAngle2AxisID] = true;
            population_.add(new Solution(lastFront.get(minAngle2AxisID)));
            remain--;
        } // for o

//		int k = 0;
//		int t = 0;
//		while (k < this.objectives_ && remain > 0) {
//			if (removed[t]== false) {
//				 population_.add(new Solution(lastFront.get(t)));
//				 removed[t] = true;
//				 k++;
//				 remain--;
//			}
//			t ++;
//		} // Add better solutions in terms of the fitness value


        /**
         * Associate each solution in the last front with a solution in the population
         */
        for (int i = 0; i < n; i++) {
            Solution solLastFront = lastFront.get(i);
            double minAng = 1.0e+30;
            int minAngID = -1;

            for (int j = 0; j < population_.size(); j++) {
                Solution solPop = population_.get(j);
                double angle = calAngle(solLastFront, solPop);
//				double angle = AngleMatrix[solLastFront.getID()][solPop.getID()];
                if (angle < minAng) {
                    minAng = angle;
                    minAngID = j;
                }
            } // for j
            angles[i] = minAng;
            index[i] = minAngID;
        } // for i

        /**
         * Niching procedure
         */
        for (int r = 0; r < remain; r++) {
            /**
             * Step 1: Find max and min angles
             */
            int maxAngleID = -1;
            double maxAngle = -1.0e+30;

            int minAglID = -1;
            double minAgl = 1.0e+30;

            for (int j = 0; j < n; j++) {
                // Find max angle
                if (removed[j] == false && angles[j] > maxAngle) {
                    maxAngle = angles[j];
                    maxAngleID = j;
                }

                // Find min angle
                if (removed[j] == false && angles[j] < minAgl) {
                    minAgl = angles[j];
                    minAglID = j;
                }
            } // for


            /**
             * Step 2: Maximum-angle-first principle
             *
             */

            if (maxAngleID != -1) {    // Not all solutions in the last front have been added
                removed[maxAngleID] = true;
                population_.add(new Solution(lastFront.get(maxAngleID)));

                // Update angles
                for (int i = 0; i < n; i++) { // For each solution in the last front
                    if (removed[i] == false) {
                        double angle = calAngle(lastFront.get(i), lastFront.get(maxAngleID));
//						 double angle = AngleMatrix[lastFront.get(i).getID()][lastFront.get(maxAngleID).getID()];
                        if (angle < angles[i]) {
                            angles[i] = angle;
                            index[i] = population_.size() - 1;
                        }
                    }//if
                }//for i

            } else {
                break;
            }

        }// for r

    }

    /**
     * @param solutionSet
     * @return
     */
    public SolutionSet readNondominatedSet(SolutionSet solutionSet) {
        SolutionSet nondominated;
        boolean[] dominated = flagsDominateNondominate(solutionSet);
        int noOfNondominated = 0;

        for (int i = 0; i < solutionSet.size(); i++) {
            if (dominated[i] == false) {
                noOfNondominated++;
            }
        }

        nondominated = new SolutionSet(noOfNondominated);

        for (int i = 0; i < solutionSet.size(); i++) {
            if (dominated[i] == false) {
                nondominated.add(solutionSet.get(i));
            }
        }

        return nondominated;
    }

    /**
     * Initialize the ideal objective vector
     *
     * @throws JMException
     * @throws ClassNotFoundException
     */
    void initIdealPoint(SolutionSet solutionSet) {
        for (int i = 0; i < objectives_; i++) {
            zp_[i] = 1.0e+30;//
        }
        for (int i = 0; i < solutionSet.size(); i++)
            updateReference(solutionSet.get(i));

    } // initIdealPoint

    /**
     * Initialize the nadir point
     *
     * @throws JMException
     * @throws ClassNotFoundException
     */
    void initNadirPoint(SolutionSet solutionSet) {
        for (int i = 0; i < objectives_; i++)
            nzp_[i] = -1.0e+30; //

        for (int i = 0; i < solutionSet.size(); i++)
            updateNadirPoint(solutionSet.get(i));

    } // initNadirPoint

    /**
     * Update the ideal objective vector
     *
     * @param indiv
     */
    void updateReference(Solution indiv) {

        for (int i = 0; i < objectives_; i++) {
            if (normalize_ == true) {
                if (indiv.getTranslatedObjectives(i) < zp_[i]) {
                    zp_[i] = indiv.getTranslatedObjectives(i);
                } // if
            } else {
                if (indiv.getObjective(i) < zp_[i]) {
                    zp_[i] = indiv.getObjective(i);
                } // if
            }
        } // for

    } // updateReference

    /**
     * Update the nadir point
     *
     * @param indiv
     */
    void updateNadirPoint(Solution indiv) {

        for (int i = 0; i < objectives_; i++) {

            if (normalize_ == true) {
                if (indiv.getTranslatedObjectives(i) > nzp_[i])
                    nzp_[i] = indiv.getTranslatedObjectives(i);
            } else {
                if (indiv.getObjective(i) > nzp_[i])
                    nzp_[i] = indiv.getObjective(i);
            }

        } // for

    } // updateNadirPoint

    /**
     * Compute the norm of a solution
     *
     * @param sol
     */
    public void computeNorm(Solution sol) {

        double norm = 0.0;
        for (int j = 0; j < objectives_; j++) {
            if (normalize_ == true) {
                norm += (sol.getTranslatedObjectives(j) - referencePoint_[j])
                        * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
            } else {
                norm += (sol.getObjective(j) - referencePoint_[j])
                        * (sol.getObjective(j) - referencePoint_[j]);
            }
        }
        norm = Math.sqrt(norm);
        sol.setDistanceToReferencePoint(norm); // This is the norm of the individual

    }// computeNorm

    /**
     * Compute the norm of a solution
     *
     * @param sol
     */
    public void computeDistance2IdealPoint(Solution sol) {

        double norm = 0.0;
        for (int j = 0; j < objectives_; j++) {
            if (normalize_ == true) {
                norm += (sol.getTranslatedObjectives(j)) * (sol.getTranslatedObjectives(j));
            } else {
                norm += (sol.getObjective(j)) * (sol.getObjective(j));
            }
        }
        norm = Math.sqrt(norm);
        sol.setDistanceToIdealPoint(norm); // This is the norm of the individual

    }// computeNorm

    /**
     * Compute the norm of a vector
     */
    public double computeNorm(double[] v) {

        double norm = 0.0;

        for (int j = 0; j < objectives_; j++) {
            norm += Math.pow(v[j] - referencePoint_[j], 2);
        } // for

        norm = Math.pow(norm, 1.0 / 2);
        return norm;
    }// computeNorm

    /**
     * Compute the norm of a vector
     */
    public double computeNormFromIdealPoint(double[] v) {

        double norm = 0.0;

        for (int j = 0; j < objectives_; j++) {
            norm += Math.pow(v[j], 2);
        } // for

        norm = Math.sqrt(norm);
        return norm;
    }// computeNormFromIdealPoint

    /**
     * Compute inner product of two solutions
     *
     * @param s1
     * @param s2
     * @return
     */
    public double ComputeInnerProduct(Solution s1, Solution s2) {

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            if (normalize_ == true) {
                innerProduct += (s1.getTranslatedObjectives(i) - referencePoint_[i])
                        * (s2.getTranslatedObjectives(i) - referencePoint_[i]);
            } else {
                innerProduct += (s1.getObjective(i) - referencePoint_[i])
                        * (s2.getObjective(i) - referencePoint_[i]);
            }
        } // for

        return innerProduct;
    } // ComputeInnerProduct

    /**
     * Compute inner product of two vectors
     *
     * @return
     */
    public double ComputeInnerProduct(double[] v1, double[] v2) {

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            innerProduct += (v1[i] - referencePoint_[i]) * (v2[i] - referencePoint_[i]);
        } // for

        return innerProduct;
    } // ComputeInnerProduct


    void computeIdealPoint(SolutionSet solutionSet) {
//		zideal = new double[objectives_];

        for (int j = 0; j < objectives_; j++) {
            zideal[j] = 1.0e+30;

            for (int i = 0; i < solutionSet.size(); i++) {
                if (solutionSet.get(i).getObjective(j) < zideal[j])
                    zideal[j] = solutionSet.get(i).getObjective(j);
            }

//			System.out.println("zideal[j] " + zideal[j] );
        }

    }

    void computeMaxPoint(SolutionSet solutionSet) {
//		zmax = new double[objectives_];

        for (int j = 0; j < objectives_; j++) {
            zmax[j] = -1.0e+30;

            for (int i = 0; i < solutionSet.size(); i++) {
                if (solutionSet.get(i).getObjective(j) > zmax[j])
                    zmax[j] = solutionSet.get(i).getObjective(j);
            }
//			System.out.println("zmax[j] " + zmax[j] );
        }
    }


    void normalizePopulationZminZmax() {
        for (int i = 0; i < union_.size(); i++) {
            Solution sol = union_.get(i);

            for (int j = 0; j < objectives_; j++) {

                double val = (sol.getObjective(j) - zideal[j])
                        / (zmax[j] - zideal[j]);

                sol.setTranslatedObjectives(val, j);
            }// for

        }

    }


    public void computeNorm(SolutionSet front) {
        // Calculate norm of each solution in front
        for (int i = 0; i < front.size(); i++) {
            Solution sol = front.get(i);

            double norm = 0.0;
            for (int j = 0; j < objectives_; j++) {
                norm += (sol.getTranslatedObjectives(j) - referencePoint_[j])
                        * (sol.getTranslatedObjectives(j) - referencePoint_[j]);
            }
            norm = Math.sqrt(norm);
            sol.setDistanceToReferencePoint(norm); // This is the norm of the individual
        }
    }

    public void computeNormFromNZP(SolutionSet front) {
        // Calculate norm of each solution in front
        for (int i = 0; i < front.size(); i++) {
            Solution sol = front.get(i);

            double norm = 0.0;
            for (int j = 0; j < objectives_; j++) {
                norm += (sol.getTranslatedObjectives(j) - nzp_[j])
                        * (sol.getTranslatedObjectives(j) - nzp_[j]);
            }
            norm = Math.sqrt(norm);
            sol.setDistanceToReferencePoint(norm); // This is the norm of the individual
        }
    }

    /**
     * Calculate the angle between two solutions
     *
     * @param s1
     * @param s2
     * @return
     */
    public double calAngle(Solution s1, Solution s2) {
        double angle = 0.0;
        double norm1 = 0.0;
        double norm2 = 0.0;

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            innerProduct += (s1.getTranslatedObjectives(i) - referencePoint_[i])
                    * (s2.getTranslatedObjectives(i) - referencePoint_[i]);
        }

        norm1 = s1.getDistanceToReferencePoint();
        norm2 = s2.getDistanceToReferencePoint();

        double value = Math.abs(innerProduct / (norm1 * norm2));
        if (value > 1.0) {
            value = 1.0;
        }

        angle = Math.acos(value);

//		Float f = new Float (angle);
//		if (f.isNaN()) {
//			System.out.println( norm1);
//		}

        return angle;
    }

    /**
     * Calculate the angle between two solutions
     *
     * @param s1
     * @param s2
     * @return
     */
    public double calAngleFromNZP(Solution s1, Solution s2) {
        double angle = 0.0;
        double norm1 = 0.0;
        double norm2 = 0.0;

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            innerProduct += (s1.getTranslatedObjectives(i) - nzp_[i])
                    * (s2.getTranslatedObjectives(i) - nzp_[i]);
        }

        norm1 = s1.getDistanceToReferencePoint();
        norm2 = s2.getDistanceToReferencePoint();

        double value = Math.abs(innerProduct / (norm1 * norm2));
        if (value > 1.0) {
            value = 1.0;
        }

        angle = Math.acos(value);
        return angle;
    }

    /**
     * Calculate the angle between a solution and a vector
     *
     * @param s1
     * @param vector
     * @return
     */
    public double calAngle(Solution s1, double[] vector) {
        double angle = 0.0;
        double norm1 = 0.0;
        double norm2 = this.computeNorm(vector);

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            innerProduct += (s1.getTranslatedObjectives(i) - referencePoint_[i])
                    * (vector[i] - referencePoint_[i]);
        }

        norm1 = s1.getDistanceToReferencePoint();

        angle = Math.acos(Math.abs(innerProduct / (norm1 * norm2)));

        return angle;
    }

    /**
     * Calculate the angle between a solution and a vector
     *
     * @param s1
     * @param vector
     * @return
     */
    public double calAngleFromIdealPoint(Solution s1, double[] vector) {
        double angle = 0.0;
        double norm1 = 0.0;
        double norm2 = computeNormFromIdealPoint(vector);

        double innerProduct = 0.0;

        for (int i = 0; i < objectives_; i++) {
            innerProduct += (s1.getTranslatedObjectives(i)) * (vector[i]);
        }

        norm1 = s1.getDistanceToIdealPoint();

        angle = Math.acos(Math.abs(innerProduct / (norm1 * norm2)));

        return angle;
    }
} // PAEAII

