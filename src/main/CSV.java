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
	public static final double[][] read(File file, char delimiter) throws NumberFormatException, IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		List<double[]> list = new ArrayList<double[]>();
		String line;
		while ((line = in.readLine()) != null) {
			String[] elements = line.split("\\s*"+delimiter+"\\s*");
			double[] row = new double[elements.length];
			for (int j = 0; j < elements.length; j ++)
				row[j] = Double.parseDouble(elements[j]);
			list.add(row);
		}
		in.close();
		return list.toArray(new double[0][]);
	}
	
}
