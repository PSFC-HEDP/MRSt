/**
 * MIT License
 *
 * Copyright (c) 2022 Justin Kunimune
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
package util;

import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;

public class InputParser {
	public String filename;
	public String implosionName;
	public int numYields;
	public int numCores;
	public IonOpticConfiguration opticsConfig;
	public DetectorConfiguration detectorConfig;
	public double uncertainty;
	public double energyBin;
	public double timeBin;
	public double tolerance;

	public InputParser(String name, String[] args) {
		// first, parse the arguments
		StringBuilder filename = new StringBuilder(name);
		this.implosionName = "og with falling temp";
		this.numYields = 1000;
		this.numCores = Math.min(10, Runtime.getRuntime().availableProcessors());
		this.opticsConfig = null;
		this.detectorConfig = null;
		this.uncertainty = 0;
		this.energyBin = 50e-3;
		this.timeBin = 20e-3;
		this.tolerance = 0.1;

		for (String arg : args) {
			if (arg.contains("=")) {
				String key = arg.substring(0, arg.indexOf('='));
				String value = arg.substring(arg.indexOf('=') + 1);
				String tagFormat = "_%.6s";
				switch (key) {
					case "implosion":
						this.implosionName = value;
						break;
					case "yields":
						this.numYields = Integer.parseInt(value);
						break;
					case "cores":
						numCores = Integer.parseInt(value);
						tagFormat = "";
						break;
					case "uncertainty":
						this.uncertainty = Double.parseDouble(value);
						tagFormat = "_%sc";
						break;
					case "energyBin":
						this.energyBin = Double.parseDouble(value)*1e-3;
						tagFormat = "_%skeV";
						break;
					case "timeBin":
						this.timeBin = Double.parseDouble(value)*1e-3;
						tagFormat = "_%sps";
						break;
					case "tolerance":
						this.tolerance = Double.parseDouble(value);
						break;
					case "width":
						if (this.detectorConfig != null)
							for (int i = 0; i < this.detectorConfig.slitWidths.length; i ++)
								detectorConfig.slitWidths[i] = Double.parseDouble(value)*1e-6;
						else
							throw new IllegalArgumentException("error! slit width was supplied before detector configuration");
						tagFormat = "_%sum";
						break;
					case "optics":
						if (value.toLowerCase().startsWith("h"))
							this.opticsConfig = IonOpticConfiguration.HIGH_EFFICIENCY;
						else if (value.toLowerCase().startsWith("m"))
							this.opticsConfig = IonOpticConfiguration.MID_EFFICIENCY;
						else if (value.toLowerCase().startsWith("l"))
							this.opticsConfig = IonOpticConfiguration.LOW_EFFICIENCY;
						else
							System.err.println("I don't know the '" + value + "' configuration");
						break;
					case "detector":
						if (value.toLowerCase().startsWith("1"))
							this.detectorConfig = DetectorConfiguration.SINGLE_STREAK_CAMERA;
						else if (value.toLowerCase().startsWith("2"))
							this.detectorConfig = DetectorConfiguration.DOUBLE_STREAK_CAMERA;
						else if (value.toLowerCase().startsWith("d"))
							this.detectorConfig = DetectorConfiguration.DOWNSCATTER_SLIT;
						else if (value.toLowerCase().startsWith("m"))
							this.detectorConfig = DetectorConfiguration.MAXIMUM_COVERAGE;
						else
							System.err.println("I don't know the '" + value + "' camera");
						break;
					default:
						System.err.println("I don't know to what '" + key + "' refers");
						break;
				}

				String tag = value.replace(" ", "");
				filename.append(String.format(tagFormat, tag));
			}
			else {
				System.err.println("I don't understand '"+arg+"'.");
			}

			this.filename = String.format("output/%s_%tF",
										  filename,
										  System.currentTimeMillis());
		}

		if (this.opticsConfig == null)
			throw new IllegalArgumentException("you need to always specify the ion optic configuration from now on.");
		if (this.detectorConfig == null)
			throw new IllegalArgumentException("you need to always specify the detector configuration from now on.");
	}
}
