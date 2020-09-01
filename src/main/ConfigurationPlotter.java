/**
 * MIT License
 * 
 * Copyright (c) 2018 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package main;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import main.CSV.COSYMapping;

/**
 * @author Justin Kunimune
 */
public class ConfigurationPlotter {
	
	private static final int N = 100;
	private static final int M = 97;
	private static final double COSY_MINIMUM_ENERGY = 10.7e6;
	private static final double COSY_MAXIMUM_ENERGY = 14.2e6;
	private static final double COSY_REFERENCE_ENERGY = 12.45e6;
	private static final double[][] PARAM_SWEEP = {
			{300e-6, 80e-6, 4e-3, 20e-3}, // og
			{150e-6, 50e-6, 2e-3, 20e-3}, // h
			{225e-6, 50e-6, 3e-3, 20e-3}, // m
			{300e-6, 80e-6, 4e-3, 20e-3}, // l
			{400e-6, 120e-6, 5e-3, 20e-3}, // x
			{300e-6, 80e-6, 3e-3, 20e-3}, // n
			{300e-6, 80e-6, 4e-3, 40e-3}, // g
			{400e-6, 60e-6, 4e-3, 20e-3}, // f
			{400e-6, 80e-6, 5e-3, 30e-3}, // b
			{400e-6, 80e-6, 4e-3, 30e-3}, // s
			{400e-6, 60e-6, 4e-3, 40e-3}, // a
			{300e-6, 20e-6, 4e-3, 20e-3},
			{300e-6, 30e-6, 4e-3, 20e-3},
			{300e-6, 40e-6, 4e-3, 20e-3},
			{300e-6, 50e-6, 4e-3, 20e-3},
			{300e-6, 60e-6, 4e-3, 20e-3},
			{300e-6, 70e-6, 4e-3, 20e-3},
			{300e-6, 90e-6, 4e-3, 20e-3},
			{300e-6, 100e-6, 4e-3, 20e-3},
			{300e-6, 110e-6, 4e-3, 20e-3},
			{300e-6, 120e-6, 4e-3, 20e-3},
			{300e-6, 80e-6, 1.0e-3, 20e-3},
			{300e-6, 80e-6, 1.5e-3, 20e-3},
			{300e-6, 80e-6, 2.0e-3, 20e-3},
			{300e-6, 80e-6, 2.5e-3, 20e-3},
			{300e-6, 80e-6, 3.0e-3, 20e-3},
			{300e-6, 80e-6, 3.5e-3, 20e-3},
			{300e-6, 80e-6, 4.5e-3, 20e-3},
			{300e-6, 80e-6, 5.0e-3, 20e-3},
			{300e-6, 80e-6, 5.5e-3, 20e-3},
			{300e-6, 80e-6, 6.0e-3, 20e-3},
			{100e-6, 80e-6, 4e-3, 20e-3},
			{150e-6, 80e-6, 4e-3, 20e-3},
			{200e-6, 80e-6, 4e-3, 20e-3},
			{250e-6, 80e-6, 4e-3, 20e-3},
			{300e-6, 80e-6, 4e-3, 20e-3},
			{350e-6, 80e-6, 4e-3, 20e-3},
			{400e-6, 80e-6, 4e-3, 20e-3},
			{450e-6, 80e-6, 4e-3, 20e-3},
	};
	
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		for (double[] params: PARAM_SWEEP) {
			MRSt mc = null;
			COSYMapping map = CSV.readCosyCoefficients(new File("data/MRSt_IRF_FP tilted_final.txt"), 3);
			double[][] cosyCoefficients = map.coefficients;
			int[][] cosyExponents = map.exponents;
			mc = new MRSt(
					Particle.D,
					3e-3,
					params[0],
					params[1],
					CSV.read(new File("data/stopping_power_deuterons.csv"), ','),
					6e0,
					params[2],
					params[3],
					COSY_MINIMUM_ENERGY,
					COSY_MAXIMUM_ENERGY,
					COSY_REFERENCE_ENERGY,
					cosyCoefficients,
					cosyExponents,
					68,
					.1,
					null); // make the simulation
			double[] E = new double[N+1];
			for (int i = 0; i <= N; i ++)
				E[i] = 13.5 + 1.*i/N;
			double[] T = new double[M+1];
			for (int i = 0; i <= M; i ++)
				T[i] = 16.2 + .1*i/M;
			double[][] S = new double[N][M];
			for (int i = 0; i < N; i ++)
				for (int j = 0; j < M; j ++)
					S[i][j] = 0;
			S[N/2][M/2] = 1;
			double[][] Sp = mc.response(E, T, S, false);
			double[] SSlicedInEnergy = new double[Sp.length];
			for (int i = 0; i < Sp.length; i ++)
				for (int j = 0; j < Sp[0].length; j ++)
					SSlicedInEnergy[i] += Sp[i][j];
			double[] SSlicedInTime = new double[Sp[0].length];
			for (int i = 0; i < Sp.length; i ++)
				for (int j = 0; j < Sp[0].length; j ++)
					SSlicedInTime[j] += Sp[i][j];
			
			double max = 0; // find time resolution
			for (int j = 0; j < Sp[0].length; j ++)
				if (SSlicedInTime[j] > max)
					max = SSlicedInTime[j];
			double left = 0, rite = 0;
			for (int j = 1; j < Sp[0].length; j ++) {
				if (SSlicedInTime[j-1] < max/2 && SSlicedInTime[j] >= max/2) {
					double c = (max/2 - SSlicedInTime[j-1])/(SSlicedInTime[j] - SSlicedInTime[j-1]);
					left = mc.getTimeAxis()[j-1]*(1-c) + mc.getTimeAxis()[j]*c;
				}
				else if (SSlicedInTime[j-1] >= max/2 && SSlicedInTime[j] < max/2) {
					double c = (max/2 - SSlicedInTime[j])/(SSlicedInTime[j-1] - SSlicedInTime[j]);
					rite = mc.getTimeAxis()[j-1]*c + mc.getTimeAxis()[j]*(1-c);
				}
			}
			double time = (rite - left)/1e-3;
			
			max = 0; // find energy resolution
			for (int i = 0; i < Sp.length; i ++)
				if (SSlicedInEnergy[i] > max)
					max = SSlicedInEnergy[i];
			left = 0;
			rite = 0;
			for (int i = 1; i < Sp.length; i ++) {
				if (SSlicedInEnergy[i-1] < max/2 && SSlicedInEnergy[i] >= max/2) {
					double c = (max/2 - SSlicedInEnergy[i-1])/(SSlicedInEnergy[i] - SSlicedInEnergy[i-1]);
					left = mc.getEnergyAxis()[i-1]*(1-c) + mc.getEnergyAxis()[i]*c;
				}
				else if (SSlicedInEnergy[i-1] >= max/2 && SSlicedInEnergy[i] < max/2) {
					double c = (max/2 - SSlicedInEnergy[i])/(SSlicedInEnergy[i-1] - SSlicedInEnergy[i]);
					rite = mc.getEnergyAxis()[i-1]*c + mc.getEnergyAxis()[i]*(1-c);
				}
			}
			double energy = (rite - left)/1e-3;
			
			System.out.printf("[%g, %g, %g, %g, %g, %g, %g, np.nan, np.nan, np.nan],\n", params[0], params[1], params[2], params[3], energy, time, mc.efficiency(14e6));
		}
	}

}
