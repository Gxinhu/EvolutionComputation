import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class python {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Process proc;
		try {
			proc = Runtime.getRuntime().exec("python /home/hu/Desktop/javausepython/src/main/java/demo.py");// 执行py文件
			//用输入输出流来截取结果
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