package jmetal.experiments.util;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

import java.util.Vector;

/**
 * Created by hxy on 2017/4/10.
 */
public class StateTools {
	public static NaturalRanking ranking = new NaturalRanking(TiesStrategy.AVERAGE);

	public static double[] friedmanRanking(double[][][] data) {
		// objective number
		int n = data.length;
		// problem number
		int k = data[0].length;
		// algorithm number
		int s = data[0][0].length;
		double[] meanRanking = new double[s];
		double[][][] r = new double[n][k][];
		for (int objective = 0; objective < n; objective++) {
			for (int problem = 0; problem < k; problem++) {
				r[objective][problem] = ranking.rank(data[objective][problem]);
			}
		}
		for (int algorithm = 0; algorithm < s; algorithm++) {
			for (int objective = 0; objective < n; objective++) {
				for (int problem = 0; problem < k; problem++) {
					meanRanking[algorithm] += r[objective][problem][algorithm];
				}
			}
			meanRanking[algorithm] /= (n * k);
		}
		return meanRanking;
	}

	public static WilcoxonSignedRankTestResult wilcoxonSignedRank(Vector xx, Vector yy, boolean isSmallerBetter) {
		Object[] x = xx.toArray();
		Object[] y = yy.toArray();

		WilcoxonSignedRankTestResult result = new WilcoxonSignedRankTestResult();
		if (x.length != y.length) {
			try {
				throw new Exception("x.length ~= y.length");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		double[] diff = new double[x.length];
		for (int i = 0; i < diff.length; i++) {
			diff[i] = (double) x[i] - (double) y[i];
		}
		double[] diffAbs = new double[x.length];
		for (int i = 0; i < diff.length; i++) {
			diffAbs[i] = Math.abs(diff[i]);
		}
		double[] ranks = ranking.rank(diffAbs);
		double Wplus = 0.0D;
		for (int i = 0; i < diff.length; ++i) {
			if (diff[i] > 0.0D) {
				Wplus += ranks[i];
			} else if (diff[i] == 0) {
				Wplus += ranks[i] / 2;
			}
		}
		int N = x.length;
		double Wminus = (double) (N * (N + 1)) / 2.0D - Wplus;
		double T = Math.min(Wminus, Wplus);
		double ES = (double) (N * (N + 1)) / 4.0D;
		double VarS = ES * ((double) (2 * N + 1) / 6.0D);
		double z_value = (T - ES - 0.5D) / Math.sqrt(VarS);
		NormalDistribution standardNormal = new NormalDistribution(0.0D, 1.0D);
		if (isSmallerBetter) {
			result.Rplus = Wminus;
			result.Rminus = Wplus;
		} else {
			result.Rplus = Wplus;
			result.Rminus = Wminus;
		}
		result.pValue = 2.0D * standardNormal.cumulativeProbability(z_value);
		return result;
	}
}

