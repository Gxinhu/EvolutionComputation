package jmetal.util.plot;

import jmetal.qualityIndicator.QualityIndicator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;

import javax.swing.*;
import java.awt.*;

/**
 * @author imssbora
 */
public class LineBeyend4d extends JFrame {

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

	public LineBeyend4d(String X, String Y, String name, double[][] datas, QualityIndicator indicator) {
		super(name);
		this.Xaxis = X;
		this.Yaxis = Y;

		this.data = datas;
		this.truePF = indicator.getTrueParetoFront();
		// Create dataset
		DefaultCategoryDataset dataset = createDataset();
		// Create chart
		JFreeChart chart = ChartFactory.createLineChart(
				name, // Chart title
				Xaxis, // X-Axis Labe
				Yaxis, // Y-Axis Label
				dataset,
				PlotOrientation.VERTICAL,
				false,                     // include legend
				false,                     // tooltips
				false
		);
		ChartPanel panel = new ChartPanel(chart);
		setContentPane(panel);
		CategoryPlot plot = chart.getCategoryPlot();
		LineAndShapeRenderer renderer = new LineAndShapeRenderer(true, false);
		for (int i = 0; i < data.length; i++) {
			renderer.setSeriesPaint(i, Color.BLUE);
		}

		plot.setRenderer(renderer);
	}

	private DefaultCategoryDataset createDataset() {

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int j = 0; j < data.length; j++) {
			for (int k = 1; k <= data[0].length; k++) {
				dataset.addValue((Number) data[j][k - 1], j + "_", k);
			}
		}
		return dataset;
	}
}