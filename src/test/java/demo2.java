import jmetal.qualityIndicator.util.MetricsUtil;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class demo2 {
	public static void main(String[] args) {

		MetricsUtil utils_ = new MetricsUtil();
		double[][] pf = utils_.readFront("/home/hu/Desktop/EvolutionComputation/10000PF/DTLZ/3d/DTLZ1.pf");
		int objective = pf[0].length;
		try {
			String[] args1 = new String[pf.length * objective + 3];
			args1[0] = "python";
			args1[1] = "/home/hu/Desktop/EvolutionComputation/src/test/java/demo2.py";
			args1[2] = String.valueOf(objective);
			int count = 3;
			for (int i = 0; i < pf.length; i++) {
				for (int j = 0; j < objective; j++) {
					args1[count] = String.valueOf(pf[i][j]);
					count++;
				}
			}
//			args1[2]="2.200047887871005001";
			Process proc = Runtime.getRuntime().exec(args1);// 执行py文件

			BufferedReader in = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			String line = null;
			while ((line = in.readLine()) != null) {
				System.out.println(line);
			}
			in.close();
			proc.waitFor();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}
