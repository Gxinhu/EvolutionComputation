//  FitnessBinaryTournament.java
//


package jmetal.operators.selection;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.DominanceComparator;

import java.util.Comparator;
import java.util.HashMap;

/**
 * This class implements an operator for binary selections using the same code
 * in FMOEA implementation
 */
public class FitnessBinaryTournament extends Selection {
  
  /**
   * dominance_ store the <code>Comparator</code> for check dominance_
   */
  private Comparator comparator_;
  
  /**
   * a_ stores a permutation of the solutions in the solutionSet used
   */
  private int a_[];
  
  /**
   *  index_ stores the actual index for selection
   */
  private int index_ = 0;
    
  /**
   * Constructor
   * Creates a new instance of the Binary tournament operator (Deb's
   * NSGA-II implementation version)
   */
  public FitnessBinaryTournament(HashMap<String, Object> parameters)
  {
  	super(parameters) ;
  	comparator_ = new DominanceComparator();              
  } // BinaryTournament2
    
  /**
  * Performs the operation
  * @param object Object representing a SolutionSet
  * @return the selected solution
  */
  public Object execute(Object object)    
  {
    SolutionSet population = (SolutionSet)object;
    if (index_ == 0) //Create the permutation
    {
      a_= (new jmetal.util.PermutationUtility()).intPermutation(population.size());
    }
            
        
    Solution solution1,solution2;
    solution1 = population.get(a_[index_]);
    index_ = (index_ + 1) % population.size();
    
    solution2 = population.get(a_[index_]);
        
    index_ = (index_ + 1) % population.size();
        
//   if (solution1.getFitness() < solution2.getFitness())
//    return solution1;
//   else if (solution2.getFitness() < solution1.getFitness())
//    return solution2;
//   else
//    if (PseudoRandom.randDouble() < 0.5)
//      return solution1;
//    else
//      return solution2;        
    
    int flag = comparator_.compare(solution1,solution2);
    if (flag == -1)
      return solution1;
    else if (flag == 1)
      return solution2;
    else if (solution1.getFitness() < solution2.getFitness())
      return solution1;
    else if (solution2.getFitness() < solution1.getFitness())
      return solution2;
    else
      if (PseudoRandom.randDouble()<0.5)
        return solution1;
      else
        return solution2;        
  } // execute
} // FitnessBinaryTournament
