/**
 * MIT License
 * <p>
 * Copyright (c) 2022 Justin Kunimune
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package util;

import physics.Analysis;
import physics.Detector.DetectorConfiguration;
import physics.IonOptics.IonOpticConfiguration;
import physics.Particle;

import java.io.File;
import java.io.IOException;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class InputParser {
	public String filename;
	public String implosionName;
	public int numRuns;
	public int numCores;
	public IonOpticConfiguration opticsConfig;
	public DetectorConfiguration detectorConfig;
	public Particle ion;
	public double shielding;
	public double uncertainty;
	public double energyBin;
	public double timeBin;
	public double offset;
	public double tolerance;
	public double yieldFactor;

	public InputParser(String name, String[] args) {
		// first, parse the arguments
		StringBuilder filename = new StringBuilder(name);
		this.implosionName = "scan base";
		this.numRuns = 1000;
		this.numCores = Math.min(10, Runtime.getRuntime().availableProcessors());
		this.opticsConfig = null;
		this.detectorConfig = null;
		this.ion = Particle.D;
		this.shielding = 293; // because of Alex Sandberg's thesis
		this.uncertainty = 0;
		this.energyBin = Analysis.DEFAULT_ENERGY_BIN;
		this.timeBin = Analysis.DEFAULT_TIME_BIN;
		this.tolerance = 0.01;
		this.yieldFactor = 1;

		for (String arg : args) {
			if (arg.contains("=")) {
				String key = arg.substring(0, arg.indexOf('='));
				String value = arg.substring(arg.indexOf('=') + 1);
				String tagFormat = "_%.8s";
				switch (key) {
					case "implosion":
						this.implosionName = value;
						break;
					case "runs":
						this.numRuns = Integer.parseInt(value);
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
					case "yield":
						this.yieldFactor = Double.parseDouble(value);
						tagFormat = "_%sx";
						break;
					case "width":
						if (this.detectorConfig != null)
							for (int i = 0; i < this.detectorConfig.slitWidths.length; i ++)
								detectorConfig.slitWidths[i] = Double.parseDouble(value)*1e-6;
						else
							throw new IllegalArgumentException("error! slit width was supplied before detector configuration");
						tagFormat = "_%sum";
						break;
					case "shielding":
						if (this.detectorConfig != null)
							this.shielding = Double.parseDouble(value);
						else
							throw new IllegalArgumentException("error! shielding was supplied before detector configuration");
						tagFormat = "_1in%s";
						break;
					case "ion":
						if (value.toLowerCase().startsWith("d"))
							this.ion = Particle.D;
						else if (value.toLowerCase().startsWith("p"))
							this.ion = Particle.P;
						else
							System.err.println("I don't know what particle '" + value + "' is.");
					case "optics":
						if (value.toLowerCase().startsWith("h"))
							this.opticsConfig = IonOpticConfiguration.HIGH_EFFICIENCY;
						else if (value.toLowerCase().startsWith("m"))
							this.opticsConfig = IonOpticConfiguration.MID_EFFICIENCY;
						else if (value.toLowerCase().startsWith("l"))
							this.opticsConfig = IonOpticConfiguration.LOW_EFFICIENCY;
						else if (value.toLowerCase().startsWith("p"))
							this.opticsConfig = IonOpticConfiguration.PERFECT;
						else
							System.err.println("I don't know the '" + value + "' configuration");
						break;
					case "detector":
						if (value.toLowerCase().startsWith("2"))
							this.detectorConfig = DetectorConfiguration.DOUBLE_STREAK_CAMERA;
						else if (value.toLowerCase().startsWith("m"))
							this.detectorConfig = DetectorConfiguration.MAXIMUM_COVERAGE;
						else if (value.toLowerCase().startsWith("j"))
							this.detectorConfig = DetectorConfiguration.JOHAN;
						else if (value.toLowerCase().startsWith("l"))
							this.detectorConfig = DetectorConfiguration.LONG_LONG;
						else if (value.toLowerCase().startsWith("d"))
							this.detectorConfig = DetectorConfiguration.DRIFT_TUBE;
						else if (value.toLowerCase().startsWith("p"))
							this.detectorConfig = DetectorConfiguration.PERFECT;
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
		}

		new File("output/").mkdirs();
		this.filename = String.format("output/%s_%tF",
		                              filename,
		                              System.currentTimeMillis());

		if (this.opticsConfig == null)
			throw new IllegalArgumentException("you need to always specify the ion optic configuration from now on.");
		if (this.detectorConfig == null)
			throw new IllegalArgumentException("you need to always specify the detector configuration from now on.");
	}


	/**
	 * idk where to put this but here is probably not correct.
	 */
	public static Logger setUpLogger(String filename) throws IOException {
		System.setProperty("java.util.logging.SimpleFormatter.format",
		                   "%1$tF %1$tT | %4$-7s | %5$s%6$s%n");
		// (I don't know why they make this so difficult)
		Formatter formatter = new SimpleFormatter();
		Logger logger = Logger.getLogger("app");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		consoleHandler.setLevel(Level.FINE);
		consoleHandler.setFormatter(formatter);
		logger.addHandler(consoleHandler);
		if (filename != null) {
			Handler logfileHandler = new FileHandler(filename + ".log");
			logfileHandler.setLevel(Level.INFO);
			logfileHandler.setFormatter(formatter);
			logger.addHandler(logfileHandler);
		}
		return logger;
	}
}
