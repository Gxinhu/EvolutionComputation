package jmetal.problems;

import java.io.*;

/**
 * 利用javacsv2.0做导入导出csv文件工具类<br/>
 *
 *
 * @author
 *
 */
public class Filetest {

	public static void main(String[] args) {
		try {
			File csv = new File("./writers.csv");//CSV文件

			BufferedWriter bw = new BufferedWriter(new FileWriter(csv, true));
			//新增一行数据
			bw.write(1);
			bw.close();
		} catch (FileNotFoundException e) {
			//捕获File对象生成时的异常
			e.printStackTrace();
		} catch (IOException e) {
			//捕获BufferedWriter对象关闭时的异常
			e.printStackTrace();
		}
	}

}
