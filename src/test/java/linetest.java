import javax.swing.*;

import org.math.plot.*;
import org.math.plot.plotObjects.BaseLabel;

import java.awt.*;

public class linetest {
	public static void main(String[] args) {

		// define your data
		double[] x = { 1, 2, 3, 4, 5, 6 };
		double[] y = { 45, 89, 6, 32, 63, 12 };

		// create your PlotPanel (you can use it as a JPanel)
		Plot2DPanel plot = new Plot2DPanel();

		// define the legend position
		plot.addLegend("SOUTH");

		// add a line plot to the PlotPanel
		plot.addScatterPlot("my plot", x, y);
		plot.setAxisLabel(0,"huxin");

		BaseLabel title = new BaseLabel("...My nice plot...", Color.RED, 0.5, 1.1);
		title.setFont(new Font("Courier", Font.BOLD, 20));
		plot.addPlotable(title);

		// put the PlotPanel in a JFrame moeads a JPanel
		JFrame frame = new JFrame("a plot panel");
		frame.setSize(600, 600);
		frame.setContentPane(plot);
		frame.setVisible(true);

	}
}