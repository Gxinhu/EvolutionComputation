/**
 * This is the code of VaEA
 * By Dr. Yi Xiang (gzhuxiang_yi@163.com)
 */
package jmetal.metaheuristics.vamaea;

import jmetal.core.*;
import jmetal.myutils.CaptureConvergence;
import jmetal.util.JMException;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;

public class VaEA extends Algorithm {

	private int populationSize_;

	private SolutionSet population_;
	SolutionSet offspringPopulation_;
	SolutionSet union_;
	int generations_;

	double alpha_;

	private double[] zp_; // 	ideal point 
	private double[] nzp_; // nadir point
	private int objectives_;

	Operator crossover_;
	Operator mutation_;
	Operator selection_;

	boolean normalize_; // Do normalization or not

	CaptureConvergence captureConvergence_; // our class 
	boolean convergenceFlag = false;

	public VaEA(Problem problem) {
		super(problem);
		objectives_ = problem.getNumberOfObjectives();
		zp_ = new double[objectives_];
		nzp_ = new double[objectives_];

		if (convergenceFlag == true) {
			captureConvergence_ = new CaptureConvergence(problem.getName(), "VaEA", problem.getNumberOfObjectives());
		}
	} // VaEA	


	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxGenerations_;

		generations_ = 0;

		maxGenerations_ = ((Integer) this.getInputParameter("maxIterations"))
				.intValue();
		populationSize_ = ((Integer) this.getInputParameter("swarmSize"))
				.intValue();
//		alpha_ = ((Double) this.getInputParameter("alpha"))
//				.doubleValue();
		alpha_ = Math.PI / 2 / (populationSize_ + 1);
		System.out.println("alpha = " + alpha_);
		normalize_ = ((Boolean) this.getInputParameter("normalize")).booleanValue();

		System.out.println("popSize = " + populationSize_);

		mutation_ = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");

		initPopulation();
		if (convergenceFlag == true) {
			captureConvergence_.runCaptureConvergence(0, population_);
		}


		while (generations_ < maxGenerations_) {
			offspringPopulation_ = new SolutionSet(populationSize_);
			Solution[] parents = new Solution[2];
			for (int i = 0; i < (populationSize_ / 2); i++) {
				if (generations_ < maxGenerations_) {
					// obtain parents

					parents = (Solution[]) selection_.execute(population_);

					Solution[] offSpring = (Solution[]) crossover_
							.execute(parents);
					mutation_.execute(offSpring[0]);
					mutation_.execute(offSpring[1]);

					problem_.evaluate(offSpring[0]);
					problem_.evaluate(offSpring[1]);


					offspringPopulation_.add(offSpring[0]);
					offspringPopulation_.add(offSpring[1]);
				} // if
			} // for	

			union_ = ((SolutionSet) population_).union(offspringPopulation_);

			population_.clear();
			new Niching(population_, union_, normalize_, populationSize_, alpha_)
					.execute();

			generations_++;

			if (convergenceFlag == true) {
				captureConvergence_.runCaptureConvergence(generations_ * populationSize_, population_);
			}

		}
		System.out.println("popSize = " + population_.size());
		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);

	}

	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);

			population_.add(newSolution);
		} // for
	} // initPopulation

	/**
	 * Initialize the ideal objective vector
	 *
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initIdealPoint() {
		for (int i = 0; i < objectives_; i++) {
			zp_[i] = 1.0e+30;
		}

		for (int i = 0; i < population_.size(); i++) {
			updateReference(population_.get(i));
		}
	} // initIdealPoint

	/**
	 * Initialize the nadir point
	 *
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initNadirPoint() {
		for (int i = 0; i < objectives_; i++) {
			nzp_[i] = -1.0e+30;
		}

		for (int i = 0; i < population_.size(); i++) {
			updateNadirPoint(population_.get(i));
		}
	} // initNadirPoint

	/**
	 * Update the ideal objective vector
	 *
	 * @param indiv
	 */
	boolean updateReference(Solution indiv) {
		boolean flag = false;

		for (int i = 0; i < objectives_; i++) {
			if (indiv.getObjective(i) < zp_[i]) {
				zp_[i] = indiv.getObjective(i);
				flag = true;
			} // if
		} // for

		return flag;
	} // updateReference

	/**
	 * Update the nadir point
	 *
	 * @param indiv
	 */
	boolean updateNadirPoint(Solution indiv) {
		boolean flag = false;

		for (int i = 0; i < objectives_; i++) {
			if (indiv.getObjective(i) > nzp_[i]) {
				nzp_[i] = indiv.getObjective(i);
			}
		}

		return flag;
	} // updateNadirPoint
} // VaEA

