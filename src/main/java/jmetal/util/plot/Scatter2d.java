package jmetal.util.plot;
/**
 *
 */

import jmetal.qualityIndicator.QualityIndicator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

import java.awt.*;
//import org.jfree.ui.Spacer;

/**
 * A simple demonstration application showing how to create a line chart using data from an
 * {@link XYDataset}.
 *
 */
public class Scatter2d extends ApplicationFrame {
	private String Xaxis;
	private String Yaxis;
	private double[][] data;
	private String name;
	private double[][] truePF;
	private boolean plotPF;
	/**
	 *
	 */
	private static final long serialVersionUID = -2318973151598624669L;

	public Scatter2d(String X, String Y, String name, double[][] datas, QualityIndicator indicator, boolean TruePF) {
		super(name);
		this.Xaxis = X;
		this.Yaxis = Y;

		this.data = datas;
		this.truePF = indicator.getTrueParetoFront().writeObjectivesToMatrix();
		this.plotPF = TruePF;

//	public Scatter2d(final String title) {

		final XYDataset dataset = createDataset();
		final JFreeChart chart = createChart(dataset);
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(800, 500));
		setContentPane(chartPanel);

	}

	/**
	 * Creates a sample dataset.
	 *
	 * @return a sample dataset.
	 */
	private XYDataset createDataset() {

		final XYSeries series1 = new XYSeries("PF");
		for (int i = 0; i < data.length; i++) {
			series1.add(data[i][0], data[i][1]);
		}
		final XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series1);
		if (plotPF) {
			final XYSeries series2 = new XYSeries("truePF");
			for (int i = 0; i < truePF.length; i++) {
				series2.add(truePF[i][0], truePF[i][1]);
			}
			dataset.addSeries(series2);
		}

		return dataset;
	}

	/**
	 * Creates a chart.
	 *
	 * @param dataset  the data for the chart.
	 *
	 * @return a chart.
	 */
	private JFreeChart createChart(final XYDataset dataset) {

		// create the chart...
		final JFreeChart chart = ChartFactory.createScatterPlot(
				name,      // chart title
				Xaxis,                      // x axis label
				Yaxis,                      // y axis label
				dataset,                  // data
				PlotOrientation.VERTICAL,
				true,                     // include legend
				true,                     // tooltips
				true                     // urls
		);

		chart.setBackgroundPaint(Color.white);
		return chart;

	}


}