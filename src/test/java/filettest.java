
import java.io.*;

public class filettest {
	public static boolean mkDirectory(String path) {
		File file = null;
		try {
			file = new File(path);
			if (!file.exists()) {
				return file.mkdirs();
			} else {
				return false;
			}
		} catch (Exception e) {
		} finally {
			file = null;
		}
		return false;
	}

	public static void main(String[] args) {
		String mkDirectoryPath = "out/hv";
		mkDirectory(mkDirectoryPath);
		try {
			File csv = new File(mkDirectoryPath + "/" + "runtimes" + ".csv");//CSV文件
			if (csv.exists()) {
				csv.delete();
				csv.createNewFile();
			} else {
				csv.createNewFile();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(csv, true));
			//新增一行数据
			bw.flush();
			bw.write("hello");
			bw.close();
			//新增一行数据
		} catch (FileNotFoundException e) {
			//捕获File对象生成时的异常
			e.printStackTrace();
		} catch (IOException e) {
			//捕获BufferedWriter对象关闭时的异常
			e.printStackTrace();
		}
		System.out.println("hello ");
	}

}
