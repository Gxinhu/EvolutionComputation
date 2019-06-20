package jmetal.metaheuristics.agmopso;

public class TestStatistics {
	private double[] inputData;

	public TestStatistics(double[] inputData) {
		this.inputData = inputData;
	}


	public double getMax() {
		if (inputData == null || inputData.length == 0) {
			return -1;
		}
		int len = inputData.length;
		double max = inputData[0];
		for (int i = 0; i < len; i++) {
			if (max < inputData[i]) {
				max = inputData[i];
			}
		}
		return max;
	}


	public double getMin() {
		if (inputData == null || inputData.length == 0) {
			return -1;
		}
		int len = inputData.length;
		double min = inputData[0];
		for (int i = 0; i < len; i++) {
			if (min > inputData[i]) {
				min = inputData[i];
			}
		}
		return min;
	}


	public double getSum() {
		if (inputData == null || inputData.length == 0) {
			return -1;
		}
		int len = inputData.length;
		double sum = 0;
		for (int i = 0; i < len; i++) {
			sum = sum + inputData[i];
		}

		return sum;

	}


	public int getCount() {
		if (inputData == null) {
			return -1;
		}

		return inputData.length;
	}


	public double getAverage() {
		if (inputData == null || inputData.length == 0) {
			return -1;
		}
		int len = inputData.length;
		double result;
		result = getSum() / len;

		return result;
	}


	public double getSquareSum() {
		if (inputData == null || inputData.length == 0) {
			return -1;
		}
		int len = inputData.length;
		double sqrsum = 0.0;
		for (int i = 0; i < len; i++) {
			sqrsum = sqrsum + inputData[i] * inputData[i];
		}


		return sqrsum;
	}


	public double getVariance() {
		int count = getCount();
		double sqrsum = getSquareSum();
		double average = getAverage();
		double result;
		result = (sqrsum - count * average * average) / count;

		return result;
	}


	public double getStandardDiviation() {
		double result;

		result = Math.sqrt(Math.abs(getVariance()));

		return result;

	}
}
