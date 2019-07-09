package jmetal.util;

import jmetal.core.Problem;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

public class createWeight {
	Problem problem_;
	int populationSize;
	double[][] lamdaVectors;

	public createWeight(Problem problem, int populationSize, double[][] lamdaVectors) {
		this.lamdaVectors = lamdaVectors;
		this.populationSize = populationSize;
		this.problem_ = problem;
	}

	public double[][] initUniformWeightwithWs() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			for (int n = 0; n < populationSize; n++) {
				double a = 1.0 * n / (populationSize - 1);
				lamdaVectors[n][0] = a;
				lamdaVectors[n][1] = 1 - a;
				nw++;
			} // for
		} // if
		else if (problem_.getNumberOfObjectives() == 3) {
			int H_ = 13;
			int i, j;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						lamdaVectors[nw][0] = (double) (1.0 * i) / H_;
						lamdaVectors[nw][1] = (double) (1.0 * j) / H_;
						lamdaVectors[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
						nw++;
					} // if
				} // for
			} // for
		} // else
		else if (problem_.getNumberOfObjectives() == 5) {
			int H_ = 6;
			int a, b, c, d;
			for (a = 0; a <= H_; a++) {
				for (b = 0; b <= H_; b++) {
					for (c = 0; c <= H_; c++) {
						for (d = 0; d <= H_; d++) {
							if (a + b + c + d <= H_) {
								lamdaVectors[nw][0] = (double) (1.0 * a) / H_;
								lamdaVectors[nw][1] = (double) (1.0 * b) / H_;
								lamdaVectors[nw][2] = (double) (1.0 * c) / H_;
								lamdaVectors[nw][3] = (double) (1.0 * d) / H_;
								lamdaVectors[nw][4] = (double) (1.0 * (H_ - a - b - c - d) / H_);
								nw++;
							}
						}
					}
				}
			}
		} else if (problem_.getNumberOfObjectives() == 8) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[120][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[36][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										if (a + b + c + d + e + f + g <= H1_) {
											lambda1[nw1][0] = (double) (1.0 * a) / H1_;
											lambda1[nw1][1] = (double) (1.0 * b) / H1_;
											lambda1[nw1][2] = (double) (1.0 * c) / H1_;
											lambda1[nw1][3] = (double) (1.0 * d) / H1_;
											lambda1[nw1][4] = (double) (1.0 * e) / H1_;
											lambda1[nw1][5] = (double) (1.0 * f) / H1_;
											lambda1[nw1][6] = (double) (1.0 * g) / H1_;
											lambda1[nw1][7] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g) / H1_);
											nw1++;
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										if (a + b + c + d + e + f + g <= H2_) {
											lambda2[nw2][0] = (double) (1.0 * a) / H2_;
											lambda2[nw2][1] = (double) (1.0 * b) / H2_;
											lambda2[nw2][2] = (double) (1.0 * c) / H2_;
											lambda2[nw2][3] = (double) (1.0 * d) / H2_;
											lambda2[nw2][4] = (double) (1.0 * e) / H2_;
											lambda2[nw2][5] = (double) (1.0 * f) / H2_;
											lambda2[nw2][6] = (double) (1.0 * g) / H2_;
											lambda2[nw2][7] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g) / H2_);
											nw2++;
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (int i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (int i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		} else if (problem_.getNumberOfObjectives() == 10) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[220][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[55][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g, h, i;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										for (h = 0; h <= H1_; h++) {
											for (i = 0; i <= H1_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H1_) {
													lambda1[nw1][0] = (double) (1.0 * a) / H1_;
													lambda1[nw1][1] = (double) (1.0 * b) / H1_;
													lambda1[nw1][2] = (double) (1.0 * c) / H1_;
													lambda1[nw1][3] = (double) (1.0 * d) / H1_;
													lambda1[nw1][4] = (double) (1.0 * e) / H1_;
													lambda1[nw1][5] = (double) (1.0 * f) / H1_;
													lambda1[nw1][6] = (double) (1.0 * g) / H1_;
													lambda1[nw1][7] = (double) (1.0 * h) / H1_;
													lambda1[nw1][8] = (double) (1.0 * i) / H1_;
													lambda1[nw1][9] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g - h - i) / H1_);
													nw1++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										for (h = 0; h <= H2_; h++) {
											for (i = 0; i <= H2_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H2_) {
													lambda1[nw2][0] = (double) (1.0 * a) / H2_;
													lambda1[nw2][1] = (double) (1.0 * b) / H2_;
													lambda1[nw2][2] = (double) (1.0 * c) / H2_;
													lambda1[nw2][3] = (double) (1.0 * d) / H2_;
													lambda1[nw2][4] = (double) (1.0 * e) / H2_;
													lambda1[nw2][5] = (double) (1.0 * f) / H2_;
													lambda1[nw2][6] = (double) (1.0 * g) / H2_;
													lambda1[nw2][7] = (double) (1.0 * h) / H2_;
													lambda1[nw2][8] = (double) (1.0 * i) / H2_;
													lambda1[nw2][9] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g - h - i) / H2_);
													nw2++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		}
//		for (int i=0;i<nw;i++){
//			for(int j=0;j<problem_.getNumberOfObjectives();j++){
//				if(lamdaVectors[i][j] == 0)
//					lamdaVectors[i][j] = 0.000001;
//			}
//		}
		if (nw != populationSize) {
			System.out.println(nw + "---" + (populationSize));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
		//applly 我也不知道的什么权值向量方法
//		RealMatrix temp = new Array2DRowRealMatrix(lamdaVectors);
//		RealVector temprow;
//		for (int i = 0; i < populationSize; i++) {
//			temprow = temp.getRowVector(i);
//			temp.setRowVector(i, temprow.mapDivide(temprow.getNorm()));
//		}
//		this.lamdaVectors = temp.getData();
		//Apply the WS-transformation on the generated weight vectors
		for (int i = 0; i < populationSize; i++) {
			double prod = 1.0, sum = 0.0;
			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
				prod = prod * lamdaVectors[i][j];
			}
			if (prod != 0.0) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / lamdaVectors[i][j];
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / lamdaVectors[i][j] / sum;
				}
			} else {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / (lamdaVectors[i][j] + 0.0000001);
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / (lamdaVectors[i][j] + 0.0000001) / sum;
				}
			}
		}
		return lamdaVectors;
	} // initUniformWeight

	public double[][] initUniformWeightwithnorm() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			for (int n = 0; n < populationSize; n++) {
				double a = 1.0 * n / (populationSize - 1);
				lamdaVectors[n][0] = a;
				lamdaVectors[n][1] = 1 - a;
				nw++;
			} // for
		} // if
		else if (problem_.getNumberOfObjectives() == 3) {
			int H_ = 13;
			int i, j;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						lamdaVectors[nw][0] = (double) (1.0 * i) / H_;
						lamdaVectors[nw][1] = (double) (1.0 * j) / H_;
						lamdaVectors[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
						nw++;
					} // if
				} // for
			} // for
		} // else
		else if (problem_.getNumberOfObjectives() == 5) {
			int H_ = 6;
			int a, b, c, d;
			for (a = 0; a <= H_; a++) {
				for (b = 0; b <= H_; b++) {
					for (c = 0; c <= H_; c++) {
						for (d = 0; d <= H_; d++) {
							if (a + b + c + d <= H_) {
								lamdaVectors[nw][0] = (double) (1.0 * a) / H_;
								lamdaVectors[nw][1] = (double) (1.0 * b) / H_;
								lamdaVectors[nw][2] = (double) (1.0 * c) / H_;
								lamdaVectors[nw][3] = (double) (1.0 * d) / H_;
								lamdaVectors[nw][4] = (double) (1.0 * (H_ - a - b - c - d) / H_);
								nw++;
							}
						}
					}
				}
			}
		} else if (problem_.getNumberOfObjectives() == 8) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[120][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[36][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										if (a + b + c + d + e + f + g <= H1_) {
											lambda1[nw1][0] = (double) (1.0 * a) / H1_;
											lambda1[nw1][1] = (double) (1.0 * b) / H1_;
											lambda1[nw1][2] = (double) (1.0 * c) / H1_;
											lambda1[nw1][3] = (double) (1.0 * d) / H1_;
											lambda1[nw1][4] = (double) (1.0 * e) / H1_;
											lambda1[nw1][5] = (double) (1.0 * f) / H1_;
											lambda1[nw1][6] = (double) (1.0 * g) / H1_;
											lambda1[nw1][7] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g) / H1_);
											nw1++;
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										if (a + b + c + d + e + f + g <= H2_) {
											lambda2[nw2][0] = (double) (1.0 * a) / H2_;
											lambda2[nw2][1] = (double) (1.0 * b) / H2_;
											lambda2[nw2][2] = (double) (1.0 * c) / H2_;
											lambda2[nw2][3] = (double) (1.0 * d) / H2_;
											lambda2[nw2][4] = (double) (1.0 * e) / H2_;
											lambda2[nw2][5] = (double) (1.0 * f) / H2_;
											lambda2[nw2][6] = (double) (1.0 * g) / H2_;
											lambda2[nw2][7] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g) / H2_);
											nw2++;
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (int i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (int i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		} else if (problem_.getNumberOfObjectives() == 10) {
			int H1_ = 3, H2_ = 2;
			int nw1 = 0, nw2 = 0;
			double[][] lambda1 = new double[220][problem_.getNumberOfObjectives()];
			double[][] lambda2 = new double[55][problem_.getNumberOfObjectives()];
			int a, b, c, d, e, f, g, h, i;
			//Generate N1
			for (a = 0; a <= H1_; a++) {
				for (b = 0; b <= H1_; b++) {
					for (c = 0; c <= H1_; c++) {
						for (d = 0; d <= H1_; d++) {
							for (e = 0; e <= H1_; e++) {
								for (f = 0; f <= H1_; f++) {
									for (g = 0; g <= H1_; g++) {
										for (h = 0; h <= H1_; h++) {
											for (i = 0; i <= H1_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H1_) {
													lambda1[nw1][0] = (double) (1.0 * a) / H1_;
													lambda1[nw1][1] = (double) (1.0 * b) / H1_;
													lambda1[nw1][2] = (double) (1.0 * c) / H1_;
													lambda1[nw1][3] = (double) (1.0 * d) / H1_;
													lambda1[nw1][4] = (double) (1.0 * e) / H1_;
													lambda1[nw1][5] = (double) (1.0 * f) / H1_;
													lambda1[nw1][6] = (double) (1.0 * g) / H1_;
													lambda1[nw1][7] = (double) (1.0 * h) / H1_;
													lambda1[nw1][8] = (double) (1.0 * i) / H1_;
													lambda1[nw1][9] = (double) (1.0 * (H1_ - a - b - c - d - e - f - g - h - i) / H1_);
													nw1++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			//Generate N2
			for (a = 0; a <= H2_; a++) {
				for (b = 0; b <= H2_; b++) {
					for (c = 0; c <= H2_; c++) {
						for (d = 0; d <= H2_; d++) {
							for (e = 0; e <= H2_; e++) {
								for (f = 0; f <= H2_; f++) {
									for (g = 0; g <= H2_; g++) {
										for (h = 0; h <= H2_; h++) {
											for (i = 0; i <= H2_; i++) {
												if (a + b + c + d + e + f + g + h + i <= H2_) {
													lambda1[nw2][0] = (double) (1.0 * a) / H2_;
													lambda1[nw2][1] = (double) (1.0 * b) / H2_;
													lambda1[nw2][2] = (double) (1.0 * c) / H2_;
													lambda1[nw2][3] = (double) (1.0 * d) / H2_;
													lambda1[nw2][4] = (double) (1.0 * e) / H2_;
													lambda1[nw2][5] = (double) (1.0 * f) / H2_;
													lambda1[nw2][6] = (double) (1.0 * g) / H2_;
													lambda1[nw2][7] = (double) (1.0 * h) / H2_;
													lambda1[nw2][8] = (double) (1.0 * i) / H2_;
													lambda1[nw2][9] = (double) (1.0 * (H2_ - a - b - c - d - e - f - g - h - i) / H2_);
													nw2++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			nw = nw1 + nw2;
			double tao = 0.5;
			for (int k = 0; k < nw2; k++) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lambda2[k][j] = (1.0 - tao) / (double) problem_.getNumberOfObjectives() + tao * lambda2[k][j];
				}
			}
			int n = 0;
			for (i = 0; i < nw1; i++) {
				lamdaVectors[n] = lambda1[i];
				n++;
			}
			for (i = 0; i < nw2; i++) {
				lamdaVectors[n] = lambda2[i];
				n++;
			}
		}
		for (int i = 0; i < nw; i++) {
			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
				if (lamdaVectors[i][j] == 0) {
					lamdaVectors[i][j] = 0.000001;
				}
			}
		}
		if (nw != populationSize) {
			System.out.println(nw + "---" + (populationSize));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
		//applly 我也不知道的什么权值向量方法
		RealMatrix temp = new Array2DRowRealMatrix(lamdaVectors);
		RealVector temprow;
		for (int i = 0; i < populationSize; i++) {
			temprow = temp.getRowVector(i);
			temp.setRowVector(i, temprow.mapDivide(temprow.getNorm()));
		}
		this.lamdaVectors = temp.getData();
		return lamdaVectors;
	} // initUniformWeight
	public double[][] initUniformWeightnorm() {
		String dataFileName;
		String dataDirectory_ = "/home/hu/Desktop/EvolutionComputation/weight/";
		dataFileName = "W" + problem_.getNumberOfObjectives() + "D_"
				+ populationSize + ".dat";

		try {
			// Open the file
			FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
					+ dataFileName);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);

			int i = 0;
			int j = 0;
			String aux = br.readLine();
			while (aux != null) {
				StringTokenizer st = new StringTokenizer(aux);
				j = 0;
				while (st.hasMoreTokens()) {
					double value = new Double(st.nextToken());
					lamdaVectors[i][j] = value;
					j++;
				}
				aux = br.readLine();
				i++;
			}
			br.close();
		} catch (Exception e) {
			System.out
					.println("initUniformWeight: failed when reading for file: "
							+ dataDirectory_ + "/" + dataFileName);
			e.printStackTrace();
		}
		RealMatrix temp = new Array2DRowRealMatrix(lamdaVectors);
		RealVector temprow;
		for (int i = 0; i < populationSize; i++) {
			temprow = temp.getRowVector(i);
			temp.setRowVector(i, temprow.mapDivide(temprow.getNorm()));
		}
		this.lamdaVectors = temp.getData();
		return lamdaVectors;
	}
	public double[][] initUniformWeightWs() {
		String dataFileName;
		String datadirectory = "/home/hu/Desktop/EvolutionComputation/weight/";
		dataFileName = "W" + problem_.getNumberOfObjectives() + "D_"
				+ populationSize + ".dat";

		try {
			// Open the file
			FileInputStream fis = new FileInputStream(datadirectory + "/"
					+ dataFileName);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);

			int i = 0;
			int j = 0;
			String aux = br.readLine();
			while (aux != null) {
				StringTokenizer st = new StringTokenizer(aux);
				j = 0;
				while (st.hasMoreTokens()) {
					double value = new Double(st.nextToken());
					lamdaVectors[i][j] = value;
					j++;
				}
				aux = br.readLine();
				i++;
			}
			br.close();
		} catch (Exception e) {
			System.out
					.println("initUniformWeight: failed when reading for file: "
							+ datadirectory + "/" + dataFileName);
			e.printStackTrace();
		}
		//Apply the WS-transformation on the generated weight vectors
		for (int i = 0; i < populationSize; i++) {
			double prod = 1.0, sum = 0.0;
			for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
				prod = prod * lamdaVectors[i][j];
			}
			if (prod != 0.0) {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / lamdaVectors[i][j];
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / lamdaVectors[i][j] / sum;
				}
			} else {
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					sum = sum + 1.0 / (lamdaVectors[i][j] + 0.0000001);
				}
				for (int j = 0; j < problem_.getNumberOfObjectives(); j++) {
					lamdaVectors[i][j] = 1.0 / (lamdaVectors[i][j] + 0.0000001) / sum;
				}
			}
		}
		return lamdaVectors;
	}
}
