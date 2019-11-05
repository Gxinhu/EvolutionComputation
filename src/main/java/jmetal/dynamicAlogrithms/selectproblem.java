package jmetal.dynamicAlogrithms;

import jmetal.core.Problem;
import jmetal.problems.dynamicProblem.DMOP.DMOP1;
import jmetal.problems.dynamicProblem.DMOP.DMOP2;
import jmetal.problems.dynamicProblem.DMOP.DMOP3;
import jmetal.problems.dynamicProblem.FDA.*;

public class selectproblem {
	Problem problem;
	int fun;
	int objectives;
	private int serverityofchange, numberofchange, t0;

	public selectproblem(Problem problem, int fun, int m, int serverityofchange, int numberofchange, int t0) {
		this.problem = problem;
		this.fun = fun;
		this.objectives = m;
		this.serverityofchange = serverityofchange;
		this.numberofchange = numberofchange;
		this.t0 = t0;
	}

	public Problem getProblem() throws ClassNotFoundException {
		switch (fun) {
			case 1: {
				problem = new FDA1("Real", 10, serverityofchange, numberofchange, t0);
				break;
			}
			case 2: {
				problem = new FDA2("Real", 13, serverityofchange, numberofchange, t0);
				break;
			}
			case 3: {
				problem = new FDA3("Real", 10, serverityofchange, numberofchange, t0);
				break;
			}
			case 4: {
				problem = new FDA4("Real", 12, serverityofchange, numberofchange, t0);
				break;
			}
			case 5: {
				problem = new FDA5("Real", 12, serverityofchange, numberofchange, t0);
				break;
			}
			case 6: {
				problem = new DMOP1("Real", 10, serverityofchange, numberofchange, t0);
				break;
			}
			case 7: {
				problem = new DMOP2("Real", 10, serverityofchange, numberofchange, t0);
				break;
			}
			case 8: {
				problem = new DMOP3("Real", 10, serverityofchange, numberofchange, t0);
				break;
			}
			default:
				throw new IllegalStateException("Unexpected value: " + fun);
		}
		return problem;
	}

}
