package jmetal.metaheuristics.heia;


import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.Ranking;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.ObjectiveComparator;

import static jmetal.metaheuristics.heia.Utils.randomPermutation;

/**
 * Implementation of NSGA-II. This implementation of NSGA-II makes use of a
 * QualityIndicator object to obtained the convergence speed of the algorithm.
 * This version is used in the paper: A.J. Nebro, J.J. Durillo, C.A. Coello
 * Coello, F. Luna, E. Alba
 * "A Study of Convergence Speed in Multi-Objective Metaheuristics." To be
 * presented in: PPSN'08. Dortmund. September 2008.
 */

public class HEIA extends Algorithm {
	/**
	 * Constructor
	 *
	 * @param problem Problem to solve
	 */
	public HEIA(Problem problem) {
		super(problem);
	} // NSGAII

	/**
	 * �����б�
	 */
	int populationSize;
	int maxEvaluations;
	int evaluations;
	int clonesize = 0;
	int num_obj = problem_.getNumberOfObjectives();
	//int[] cro_type = new int[3];
	//int[][] neighborsets=new int[clonesize][neighborsize];

	QualityIndicator indicators; // QualityIndicator object
	//int requiredEvaluations; // Use in the example of use of the
	// indicators object (see below)
	double Cro_pro = 0.5;
	Ranking ranking;
	SolutionSet population;
	SolutionSet offspringPopulation;
	SolutionSet union;
	SolutionSet clonepopulation;
	SolutionSet Archive;
	SolutionSet front = null;
	SolutionSet SBXPop, DEPop;

	Operator mutationOperator;
	Operator crossoverOperator, crossoverOperator2;
	//Operator selectionOperator;
	Operator cloneoperator;

	Distance distance = new Distance();


	/**
	 * Runs the HEIA algorithm.
	 *
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 * solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {

		// Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize"))
				.intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations"))
				.intValue();
		indicators = (QualityIndicator) getInputParameter("indicators");
		clonesize = populationSize / 5;

		// Initialize the variables
		clonepopulation = new SolutionSet(clonesize);
		population = new SolutionSet(populationSize);
		Archive = new SolutionSet(populationSize);
		offspringPopulation = new SolutionSet(populationSize);
		SBXPop = new SolutionSet(populationSize);
		DEPop = new SolutionSet(populationSize);

		//int[][] neighborhood_= new int[populationSize][20];
		evaluations = 0;
		//requiredEvaluations = 0;

		// Read the operators
		cloneoperator = operators_.get("clone");
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover"); //DE
		crossoverOperator2 = operators_.get("crossover2");//SBX
		//selectionOperator = operators_.get("selection");

		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			evaluations++;
			population.add(newSolution);
		} // for

		ranking = new Ranking(population);
		Distance distance_ = new Distance();
		front = ranking.getSubfront(0);
		front.sort(new ObjectiveComparator(0));        // ���յ�һ��Ŀ�������
		for (int k = 0; k < front.size(); k++) {    // archive�б����֧��⣬�ھӰ��յ�һ��Ŀ�������Ƴ̶ȼ��㡣
			front.get(k).setclone_num(k);
			Archive.add(front.get(k));
		} // for
		distance.crowdingDistanceAssignment(front,
				problem_.getNumberOfObjectives());
		front.sort(new CrowdingComparator());
		for (int k = 0; k < clonesize && k < front.size(); k++) {
			clonepopulation.add(front.get(k));
		} // for

		//Archive = front;
		//selectforclone(population,clonepopulation,clonesize);

		// Generations
		while (evaluations < maxEvaluations) {
			offspringPopulation = (SolutionSet) cloneoperator.execute(clonepopulation);
			SBXPop.clear();
			DEPop.clear();
			for (int i = 0; i < offspringPopulation.size(); i++) {
				if (Cro_pro < PseudoRandom.randDouble()) {
					SBXPop.add(offspringPopulation.get(i));
				} else {
					DEPop.add(offspringPopulation.get(i));
				}
			}
			DEUpdate();
			SBXUpdate();
			ArchiveUpdate();

		} // while

		return Archive;
	} // execute


	// Archive���¿�¡
	public void ArchiveUpdate() {
		// Create the solutionSet union of solutionSet and offSpring
		union = Archive.union(DEPop);
		union = union.union(SBXPop);
		union.Suppress();
		//population=offspringPopulation;
		// Ranking the union
		ranking = new Ranking(union);

		Archive.clear();
		clonepopulation.clear();

		front = ranking.getSubfront(0);

		distance.crowdingDistanceAssignment(front,
				problem_.getNumberOfObjectives());

		front.Suppress();
		// Remain is less than front(index).size, insert only the best one
		front.sort(new CrowdingComparator());

		while (front.size() > populationSize) {
			front.remove(front.size() - 1);
			distance.crowdingDistanceAssignment(front,
					problem_.getNumberOfObjectives());
			front.sort(new CrowdingComparator());
		}
		front.sort(new ObjectiveComparator((int) (PseudoRandom.randDouble() * num_obj)));
		for (int k = 0; k < front.size(); k++) {
			front.get(k).setclone_num(k);
			Archive.add(front.get(k));
		} // for
		front.sort(new CrowdingComparator());
		for (int k = 0; k < clonesize && k < front.size(); k++) {
			clonepopulation.add(front.get(k));
		} // for
	}


	// DEPop ����
	public void DEUpdate() throws JMException {
		Solution[] offSpring = new Solution[2];
		for (int i = 0; i < DEPop.size(); i++) {
			if (evaluations < maxEvaluations) {
				Solution[] parents = new Solution[3];
				parents[2] = DEPop.get(i);
				if (clonepopulation.size() < 20) {
					if (clonepopulation.size() > 1) {
						/*int seleted=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* clonepopulation.size());
						
						int seleted2=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* (clonepopulation.size()-1));
						if(seleted2>=seleted)
							seleted2++;*/
						int[] permutation = new int[clonepopulation.size()];
						randomPermutation(permutation, clonepopulation.size());
						int seleted = permutation[0];
						int seleted2 = permutation[1];
						parents[0] = clonepopulation.get(seleted);
						parents[1] = clonepopulation.get(seleted2);
					} else {
						parents[0] = clonepopulation.get(0);
						parents[1] = clonepopulation.get(0);
					}
				} else {
					if (0.1 < PseudoRandom.randDouble()) {        // 0.9�ĸ���
						int neighbors = DEPop.get(i).getclone_num();
						/*int seleted=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* 20);
						int seleted2=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* 19);
						if(seleted2>=seleted)
							seleted2++;*/
						int[] permutation = new int[20];
						randomPermutation(permutation, 20);
						int seleted = permutation[0];
						int seleted2 = permutation[1];
						if (neighbors < 10) {
							parents[1] = Archive.get(seleted);
							parents[0] = Archive.get(seleted2);
						} else if (neighbors > (Archive.size() - 10)) {
							parents[1] = Archive.get(Archive.size() - 20 + seleted);
							parents[0] = Archive.get(Archive.size() - 20 + seleted2);
						} else {
							parents[1] = Archive.get(neighbors - 10 + seleted);
							parents[0] = Archive.get(neighbors - 10 + seleted2);
						}
						//parents[1] = Archive.get(neighborhood_[neighbors][seleted]);
						//parents[0] = Archive.get(neighborhood_[neighbors][seleted2]);
						//int[] permutation = new int[20];
						//Utils.randomPermutation(permutation, 20);
					} else {
						/*int seleted=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* clonepopulation.size());
						
						int seleted2=(int) Math
    							.floor(PseudoRandom.randDouble()
    									* (clonepopulation.size()-1));
						if(seleted2>=seleted)
							seleted2++;*/
						int[] permutation = new int[clonepopulation.size()];
						randomPermutation(permutation, clonepopulation.size());
						int seleted = permutation[0];
						int seleted2 = permutation[1];
						parents[0] = clonepopulation.get(seleted);
						parents[1] = clonepopulation.get(seleted2);
					}
				}
				offSpring[0] = (Solution) crossoverOperator.execute(new Object[]{
						parents[2], parents});
				mutationOperator.execute(offSpring[0]);
				problem_.evaluate(offSpring[0]);
				problem_.evaluateConstraints(offSpring[0]);
				DEPop.replace(i, offSpring[0]);
				evaluations += 1;
			}
		}
	}

	//SBXPop ����
	public void SBXUpdate() throws JMException {
		Solution[] offSpring = new Solution[2];
		for (int i = 0; i < SBXPop.size(); i++) {
			if (evaluations < maxEvaluations) {
				Solution[] parents = new Solution[2];
				parents[0] = SBXPop.get(i);
				parents[1] = clonepopulation.get((int) Math
						.floor(PseudoRandom.randDouble()
								* clonepopulation.size()));
				offSpring = (Solution[]) crossoverOperator2
						.execute(parents);
				mutationOperator.execute(offSpring[0]);
				problem_.evaluate(offSpring[0]);
				problem_.evaluateConstraints(offSpring[0]);
				SBXPop.replace(i, offSpring[0]);
				evaluations += 1;
			}
		}
	}


} 
