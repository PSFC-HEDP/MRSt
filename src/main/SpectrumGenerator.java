/**
 * MIT License
 * 
 * Copyright (c) 2020 Justin Kunimune
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

/**
 * @author Justin Kunimune
 *
 */
public class SpectrumGenerator {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {
		double[][] thing;
		double[] eBins;
		thing = CSV.read(new File("data/Yn-rR-Ti_150327_16p26 - Yn-rR-Ti_150327_16p26.csv"), ',', 1);
		eBins = CSV.readColumn(new File("data/Energy bins.txt"));
		
		double[] time = new double[thing.length];
		double[] ρR = new double[thing.length];
		double[] Yn = new double[thing.length];
		double[] Ti = new double[thing.length];
		double[] zero = new double[thing.length];
		for (int i = 0; i < thing.length; i ++) {
			time[i] = thing[i][0];
			Yn[i] = thing[i][1]*1e3;
			Ti[i] = thing[i][2];
			ρR[i] = thing[i][4] + thing[i][5];
			zero[i] = 0;
		}
		double[] tBins = new double[time.length + 1];
		tBins[0] = (3*time[0] - time[1])/2.;
		for (int i = 1; i < time.length; i ++)
			tBins[i] = (time[i-1] + time[i])/2.;
		tBins[time.length] = (3*time[time.length-1] - time[time.length-2])/2.;
		double[][] spectrum = MRSt.generateSpectrum(Yn, Ti, zero, zero, ρR, zero, eBins, tBins, null);
		
		CSV.writeColumn(tBins, new File("data/Time bins.txt"));
		CSV.write(spectrum, new File("data/spectrum.txt"), '\t');
	}

}
