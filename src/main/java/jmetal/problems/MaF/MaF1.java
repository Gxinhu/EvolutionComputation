//  MaF1.java 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.problems.MaF;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;

/**
 * Class representing problem MaF1
 */
public class MaF1 extends Problem {
	/**
	 * Creates a default MaF1 problem (12 variables and 3 objectives)
	 *
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF1(String solutionType) throws ClassNotFoundException {
		this(solutionType, 12, 3);
	} // DTLZ1

	/**
	 * Creates a MaF1 problem instance
	 *
	 * @param numberOfVariables
	 *            Number of  
	 * @param numberOfObjectives
	 *            Number of objective functions
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public MaF1(String solutionType, Integer numberOfVariables,
	            Integer numberOfObjectives) {
		numberOfVariables_ = numberOfVariables;
		numberOfObjectives_ = numberOfObjectives;
		numberOfConstraints_ = 0;
		problemName_ = "MaF1";

		lowerLimit_ = new double[numberOfVariables_];
		upperLimit_ = new double[numberOfVariables_];
		for (int var = 0; var < numberOfVariables; var++) {
			lowerLimit_[var] = 0.0;
			upperLimit_[var] = 1.0;
		} // for

		if (solutionType.compareTo("BinaryReal") == 0)
			solutionType_ = new BinaryRealSolutionType(this);
		else if (solutionType.compareTo("Real") == 0)
			solutionType_ = new RealSolutionType(this);
		else {
			System.out.println("Error: solution type " + solutionType
					+ " invalid");
			System.exit(-1);
		}
	}

	/**
	 * Evaluates a solution
	 *
	 * @param solution
	 *            The solution to evaluate
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		Variable[] gen = solution.getDecisionVariables();

		double[] x = new double[numberOfVariables_];
		double[] f = new double[numberOfObjectives_];
		int k = numberOfVariables_ - numberOfObjectives_ + 1;

		for (int i = 0; i < numberOfVariables_; i++)
			x[i] = gen[i].getValue();

		double g = 0.0;
		for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
			g += (x[i] - 0.5) * (x[i] - 0.5);

		for (int i = 0; i < numberOfObjectives_; i++)
			f[i] = (1.0 + g);

		for (int i = 0; i < numberOfObjectives_; i++) {
			double value = 1.0;
			for (int j = 0; j < numberOfObjectives_ - (i + 1); j++) {
				value *= x[j];
			}
			if (i != 0) {
				int aux = numberOfObjectives_ - (i + 1);
				value *= 1 - x[aux];
			} // if
			value = 1.0 - value;
			f[i] *= value;
		}// for

		for (int i = 0; i < numberOfObjectives_; i++) {
			solution.setObjective(i, f[i]);
		}
	} // evaluate

}//MaF1
