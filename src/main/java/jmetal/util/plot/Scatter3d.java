package jmetal.util.plot;

import jmetal.qualityIndicator.QualityIndicator;
import org.math.plot.Plot3DPanel;
import org.math.plot.plotObjects.BaseLabel;

import javax.swing.*;
import java.awt.*;

public class Scatter3d {
	private String Xaxis;
	private String Yaxis;
	private double[][] data;
	private String name;
	double[][] truePF;
	boolean plotPF;

	public Scatter3d(String X, String Y, String name, double[][] datas, QualityIndicator indicator, boolean TurePF) {
		this.Xaxis = X;
		this.Yaxis = Y;
		this.name = name;
		this.data = datas;
		this.truePF = indicator.getTrueParetoFront().writeObjectivesToMatrix();
		this.plotPF = TurePF;
	}

	public void plot() {
		double[] x = new double[data.length];
		double[] y = new double[data.length];
		double[] z = new double[data.length];
		for (int i = 0; i < data.length; i++) {
			x[i] = data[i][0];
			y[i] = data[i][1];
			z[i] = data[i][2];
		}
		// create your PlotPanel (you can use it as a JPanel)
		Plot3DPanel plot = new Plot3DPanel();

		// legend at SOUTH
		plot.addLegend("SOUTH");

		// add the histogram (50 slices) of x to the PlotPanel
		if (plotPF == true) {
			plot.addScatterPlot("TruePF", truePF);

		}
		plot.addScatterPlot("getPF", x, y, z);
		// add a title
		BaseLabel title = new BaseLabel(name, Color.RED, 0.5, 1.1, 1.1);
		title.setFont(new Font("Courier", Font.BOLD, 20));
		plot.addPlotable(title);

		// change name of axes
//		plot.setAxesLabels("<X>", "frequency");
		plot.setAxisLabel(0, this.Xaxis);
		plot.setAxisLabel(1, this.Yaxis);

		// put the PlotPanel in a JFrame moeads a JPanel
		JFrame frame = new JFrame("a plot panel");
		frame.setSize(1000, 1000);
		frame.setContentPane(plot);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	}
}
