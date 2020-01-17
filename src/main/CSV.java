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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * a class for reading the weird input text files that go with this.
 * 
 * @author Justin Kunimune
 */
public class CSV {
	
	/**
	 * read a simple CSV file, with any standard line break character, and return its contents
	 * as a double matrix. file must end with a line break, and elements must be parsable as
	 * doubles. whitespace adjacent to delimiters will be stripped.
	 * @param file the CSV file to open
	 * @param delimiter the delimiting character, usually ',', sometimes '\t', occasionally '|'
	 * @return I think the return value is pretty self-explanatory.
	 * @throws IOException if file cannot be found or permission is denied
	 * @throws NumberFormatException if elements are not parsable as doubles
	 */
	public static final double[][] read(File file, char delimiter)
			throws NumberFormatException, IOException {
		BufferedReader in = null;
		List<double[]> list;
		try {
			in = new BufferedReader(new FileReader(file));
			list = new ArrayList<double[]>();
			String line;
			while ((line = in.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty())
					break;
				String[] elements = line.split("\\s*"+delimiter+"\\s*");
				double[] row = new double[elements.length];
				for (int j = 0; j < elements.length; j ++)
					row[j] = Double.parseDouble(elements[j]);
				list.add(row);
			}
		} finally {
			try {
				in.close();
			} catch (NullPointerException e) {}
		}
		return list.toArray(new double[0][]);
	}
	
	/**
	 * read a CSV file where there is only one column, bypassing the need for a multi-
	 * dimensional array. values will be separated by line breaks and adjacent whitespace
	 * alone, and will be returned in a 1D array.
	 * @param file the CSV file to open
	 * @return 1D array of values from the list
	 * @throws IOException if file cannot be found or permission is denied
	 * @throws NumberFormatException if elements are not parsable as doubles
	 */
	public static final double[] readColumn(File file)
			throws NumberFormatException, IOException {
		return readColumn(file, '\n', 0);
	}
	
	/**
	 * read a CSV file and return a single column as a 1D array.
	 * @param file the CSV file to open
	 * @param delimiter the delimiting character, usually ',', sometimes '\t', occasionally '|'
	 * @param j index of the desired column
	 * @return 1D array of values in column j
	 * @throws IOException if file cannot be found or permission is denied
	 * @throws NumberFormatException if elements are not parsable as doubles
	 */
	public static final double[] readColumn(File file, char delimiter, int j)
			throws NumberFormatException, IOException {
		double[][] table = read(file, delimiter);
		double[] out = new double[table.length];
		for (int i = 0; i < table.length; i ++)
			out[i] = table[i][j];
		return out;
	}
	
	/**
	 * read a COSY-generated file and return the coefficients as a double matrix.
	 * @param file the COSY file to open
	 * @return yep
	 * @throws IOException if file cannot be found or permission is denied
	 * @throws NumberFormatException if elements are not parsable as doubles
	 */
	public static final double[][] readCosyCoefficients(File file)
			throws NumberFormatException, IOException {
		double[][] fullFile = read(file, ' '); // read it normally
		double[][] coefs = new double[fullFile.length][fullFile[0].length - 1];
		for (int i = 0; i < fullFile.length; i ++) {
			System.arraycopy(fullFile[i], 0, coefs[i], 0, coefs[i].length); // but then remove the last column
		}
		return coefs;
	}

	/**
	 * read a COSY-generated file and return the exponents as an int matrix.
	 * @param file the COSY file to open
	 * @return yep
	 * @throws IOException if file cannot be found or permission is denied
	 * @throws NumberFormatException if elements are not parsable as doubles
	 */
	public static final int[][] readCosyExponents(File file)
			throws NumberFormatException, IOException {
		double[][] fullFile = read(file, ' '); // read it normally
		int[][] exps = new int[fullFile.length][9];
		for (int i = 0; i < fullFile.length; i ++) {
			int code = (int) fullFile[i][fullFile[i].length-1]; // but then we only care about the last column
			for (int j = exps[i].length-1; j >= 0; j --) {
				exps[i][j] = code%10; // take its digits one at a time
				code /= 10; // each digit is a column
			} // I know this is inefficient since it parsed the column as a double just to then cast it to an int and extract its digits, but this was what I came up with to minimize duplicated code.
		}
		return exps;
	}
	
}
