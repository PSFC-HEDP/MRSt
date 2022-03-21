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
package app;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.logging.StreamHandler;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import physics.Analysis;
import physics.Detector.DetectorConfiguration;
import physics.Particle;
import physics.SpectrumGenerator;
import util.COSYMapping;
import util.CSV;
import physics.Analysis.ErrorMode;
import util.Math2;
import util.PythonPlot;
import util.Spinner;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class SpectrumViewer extends Application {
	
	private static final Particle ION = Particle.D;
	private static final double CENTRAL_E = 12.45;
	private static final int SPACING_0 = 16;
	private static final int SPACING_1 = 10;
	private static final int SPACING_2 = 4;

	private static final String DEFAULT_TRAJECTORY = "og with falling temp";
	
	private Spinner<Double> foilDistance;
	private Spinner<Double> foilWidth;
	private Spinner<Double> foilHeight;
	private Spinner<Double> foilThickness;
	private Spinner<Double> apertureDistance;
	private Spinner<Double> apertureWidth;
	private Spinner<Double> apertureHeight;
	private Spinner<Double> focalPlaneTilt;
	private ChoiceBox<Integer> order;
	private Spinner<Double> yieldFactor;
	private CheckBox errorBars;
	private CheckBox reuseMatrix;
	
	private COSYMapping cosyMapping;
	private double[] timeBins;
	private double[] energyBins;
	private double[][] spectrum;
	private String spectrumName;
	private Logger logger;
	
	
	/**
	 * build the GUI and display it.
	 */
	public void start(Stage stage) throws NumberFormatException {
		GridPane leftPane = new GridPane();
		leftPane.setHgap(SPACING_1);
		leftPane.setVgap(SPACING_1);
		int row = 0;
		
		this.foilDistance = new Spinner<Double>(0.5, 10.0, 3.0, 0.1);
		leftPane.add(new Label("Foil distance"), 0, row);
		leftPane.add(foilDistance, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilWidth = new Spinner<Double>(0.1, 2.0, 0.8, 0.01); // TODO maybe throw a warning if the radius >~ the distance
		leftPane.add(new Label("Foil width"), 0, row);
		leftPane.add(foilWidth, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilHeight = new Spinner<Double>(0.1, 2.0, 0.8, 0.01); // TODO maybe throw a warning if the radius >~ the distance
		leftPane.add(new Label("Foil height"), 0, row);
		leftPane.add(foilHeight, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilThickness = new Spinner<Double>(5., 500., 90., 5.);
		leftPane.add(new Label("Foil thickness"), 0, row);
		leftPane.add(foilThickness, 1, row);
		leftPane.add(new Label("μm"), 2, row);
		row ++;
		
		this.apertureDistance = new Spinner<Double>(0.01, 10.0, 6.0, 0.1);
		leftPane.add(new Label("Aper. distance"), 0, row);
		leftPane.add(apertureDistance, 1, row);
		leftPane.add(new Label("m"), 2, row);
		row ++;
		
		this.apertureWidth = new Spinner<Double>(0.1, 50.0, 5.0, 1.0);
		leftPane.add(new Label("Aper. width"), 0, row);
		leftPane.add(apertureWidth, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.apertureHeight = new Spinner<Double>(0.1, 50.0, 20.0, 1.0);
		leftPane.add(new Label("Aper. height"), 0, row);
		leftPane.add(apertureHeight, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.focalPlaneTilt = new Spinner<Double>(0.0, 89.9, 66.586, 5.0);
		leftPane.add(new Label("F. plane angle"), 0, row);
		leftPane.add(focalPlaneTilt, 1, row); // XXX: this is actually effectively a streak camera configuration dropdown
		leftPane.add(new Label("°"), 2, row);
		row ++;
		
		this.order = new ChoiceBox<Integer>(FXCollections.observableArrayList(1, 2, 3));
		order.setValue(3); // TODO when this changes it should reload the cosy map
		leftPane.add(new Label("Order"), 0, row);
		leftPane.add(order, 1, row);
		row ++;
		
		leftPane.add(chooseFileWidget("COSY map file:", stage, "MRSt_IRF_FP_70deg.txt",
				(file) -> {
					this.cosyMapping = CSV.readCosyCoefficients(file, order.getValue());
					this.cosyMapping.setConfig(ION, CENTRAL_E); // someday there may be COSY matrices that are bilt for other particles, or other energies.  that day is not today.
				}), 0, row, 3, 1);

		VBox rightPane = new VBox(SPACING_1);
		
		rightPane.getChildren().add(chooseFileWidget("Energy bin file:", stage, "energy.txt",
				(file) -> {
					this.energyBins = CSV.readColumn(file);
				}));
		
		rightPane.getChildren().add(chooseFileWidget("Time bin file:", stage, "time "+DEFAULT_TRAJECTORY+".txt",
				(file) -> {
					this.timeBins = CSV.readColumn(file);
				}));
		
		rightPane.getChildren().add(chooseFileWidget("Spectrum file:", stage, "spectrum "+DEFAULT_TRAJECTORY+".txt",
				(file) -> {
					this.spectrum = CSV.read(file, '\t');
					if (file.getName().startsWith("spectrum "))
						this.spectrumName = file.getName().substring(9, file.getName().length()-4);
					else
						this.spectrumName = "-";
				}));

		this.yieldFactor = new Spinner<Double>(1e-2, 1e+6, 100, 10);
		yieldFactor.setEditable(true);
		GridPane container = new GridPane();
		container.setHgap(SPACING_1);
		container.add(new Label("Yield factor"), 0, 0);
		container.add(yieldFactor, 1, 0);
		container.add(new Label("%"), 2, 0);
		rightPane.getChildren().add(container);

		this.errorBars = new CheckBox("Compute error bars");
		this.errorBars.setSelected(true);
		rightPane.getChildren().add(errorBars);

		this.reuseMatrix = new CheckBox("Reuse last transfer matrix");
		this.reuseMatrix.setSelected(false);
		rightPane.getChildren().add(reuseMatrix);

		Button execute = new Button("Compute!");
		execute.setOnAction((event) -> {
			if (cosyMapping == null)
				logger.severe("Please select a COSY map file.");
			else if (energyBins == null)
				logger.severe("Please select an energy bin file.");
			else if (timeBins == null)
				logger.severe("Please select a time bin file.");
			else if (spectrum == null)
				logger.severe("Come on, man. You're so close. You need a spectrum to go with those bins.");
			else {
				new Thread(() -> {
					double[] eBins, tBins;
					double[][] spec;
					try {
						eBins = energyBins.clone(); // save the current values of these spectra
						tBins = timeBins.clone();
						if (spectrum.length != eBins.length-1 || spectrum[0].length != tBins.length-1) {
							logger.info("interpreting weird spectrum file");
							spec = SpectrumGenerator.interpretSpectrumFile(tBins, eBins, spectrum); // deal with the necessary differentiation etc
						}
						else
							spec = deepClone(spectrum);
						SpectrumGenerator.modifySpectrum(spec, yieldFactor.getValue()/100.*Math2.sum(spec));
					} catch (ArrayIndexOutOfBoundsException e) {
						logger.severe("Invalid input spectrum file.");
						return;
					}
					
					logger.log(Level.INFO, "running fit on spectrum with yield factor = "+yieldFactor.getValue()/100);
					Analysis mc;
					try {
						mc = new Analysis(
								foilDistance.getValue()*1e-3,
								foilWidth.getValue()*1e-3,
								foilHeight.getValue()*1e-3,
								foilThickness.getValue()*1e-6,
								apertureDistance.getValue()*1e0,
								apertureWidth.getValue()*1e-3,
								apertureHeight.getValue()*1e-3,
								cosyMapping,
								(focalPlaneTilt.getValue() == 0) ?
									DetectorConfiguration.SINGLE_STREAK_CAMERA :
							  		DetectorConfiguration.DOWNSCATTER_SLIT,
								0,
								reuseMatrix.isSelected(),
								logger); // make the simulation

						double dispersion = mc.computeDispersion();
						double skew = mc.computeTimeSkew();
						logger.info(String.format("Dispersion: %.2f keV/mm", dispersion));
						logger.info(String.format("Time skew:  %.2f ps/keV", skew));
						logger.info(String.format("Efficiency: %.4g", mc.efficiency()));
//						double[] res = mc.computeResolution(14.);
//						logger.info(String.format("Energy res: %.2f keV", res[0]));
//						logger.info(String.format("Time res:   %.2f ps", res[1]));

						mc.respondAndAnalyze(
								eBins,
								tBins,
								spec,
								errorBars.isSelected() ? ErrorMode.HESSIAN : ErrorMode.STATISTICS); // and run it!
						
					} catch (Exception e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
						return;
					}
					
					double[][] smallSpec = Math2.downsample(tBins, eBins, spec, mc.getTimeBins(), mc.getEnergyBins());
					try { // send the data to python for plotting
						PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), smallSpec,
						                       "Original neutron spectrum");
						PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getCorrectedSpectrum(),
						            "Synthetic deuteron spectrum");
						PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getInferredSpectrum(),
						            "Fitted neutron spectrum");
						PythonPlot.plotHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getFittedSpectrum(),
						            "Fitted deuteron spectrum");
						PythonPlot.plotLines(spectrumName,
								mc.getTimeAxis(), "Time (ns)",
						 		mc.getNeutronYield(), mc.getNeutronYieldError(), "Yn (10^15/ns)",
								mc.getIonTemperature(), mc.getIonTemperatureError(), "Ti (keV)",
								mc.getArealDensity(), mc.getArealDensityError(), "ρR (g/cm^2)"
						);
						PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getEnergyBins(), mc.getCorrectedSpectrum(), mc.getFittedSpectrum(),
								"Time", "Energy (MeV)", "Synthetic deuteron spectrum", "Fitted deuteron spectrum");
						PythonPlot.compareHeatmap(mc.getTimeBins(), mc.getEnergyBins(), smallSpec, mc.getInferredSpectrum(),
								"Time", "Energy (MeV)", "Original neutron spectrum", "Fitted neutron spectrum");
					} catch (IOException e) {
						logger.log(Level.SEVERE, "Could not access plotting scripts and/or plots", e);
					}
				}).start();
			}
		});
		rightPane.getChildren().add(execute);
		
		final TextArea console = new TextArea();
		console.setEditable(false);
		console.setPrefWidth(400);
		console.setFont(Font.font("Monospace"));
		rightPane.getChildren().add(console);

		System.setProperty("java.util.logging.SimpleFormatter.format",
				"%1$tF %1$tT | %4$-7s | %5$s%6$s%n");
		logger = Logger.getLogger("app");
		logger.setUseParentHandlers(false);
		logger.setLevel(Level.ALL);
		StreamHandler guiHandler = new StreamHandler() {
			public void publish(LogRecord record) {
				Platform.runLater(() -> {
					console.appendText(String.format("%7s: %s\n",
							record.getLevel().toString(), record.getMessage()));
				});
			}
		};
		guiHandler.setLevel(Level.FINE);
		logger.addHandler(guiHandler);
		ConsoleHandler commandlineHandler = new ConsoleHandler();
		commandlineHandler.setLevel(Level.ALL);
		logger.addHandler(commandlineHandler);
		
		StackPane root = new StackPane();
		root.setPadding(new Insets(SPACING_0));
		root.getChildren().add(new HBox(SPACING_0, leftPane, rightPane));
		Scene scene = new Scene(root);
		
		stage.setTitle("Spectrum viewer");
		stage.setScene(scene);
		stage.show();
	}
	
	
	/**
	 * create a consistent-looking file selection thingy with a button and label, and bind it
	 * to the given File.
	 * @param title the title of the window
	 * @param initialFilename the filename of the file that it will try to load initially
	 * @param action the action to take with the file once it's chosen
	 * @return the button and label of the file picker
	 */
	private static Region chooseFileWidget(String title, Stage stage, String initialFilename,
			Callback action) {
		Label label = new Label();
		Button button = new Button("Chose file…");
		FileChooser fileChooser = new FileChooser();
		fileChooser.setTitle("Load "+title);
		fileChooser.setInitialDirectory(new File("input/"));
		fileChooser.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("Data files", "*.csv", "*.tsv", "*.txt"),
				new FileChooser.ExtensionFilter("All files", "*.*"));
		button.setOnAction((event) -> {
			final File chosen = fileChooser.showOpenDialog(stage);
			if (chosen != null) {
				label.setText(chosen.getName());
				try {
					action.process(chosen);
				} catch (IOException e) {
					System.err.println("Could not open "+chosen.getName()); // TODO remind myself how Alerts work
				} catch (Exception e) {
					System.err.println("There was a problem opening "+chosen.getName());
					e.printStackTrace();
				}
			}
		});
		HBox output = new HBox(SPACING_2, button, label);
		label.setMaxWidth(250 - button.getPrefWidth());
		output.setPrefWidth(250);
		output.setAlignment(Pos.CENTER_LEFT);
		label.setText(initialFilename);
		try {
			action.process(new File("input/"+initialFilename));
		} catch (IOException e) {
			label.setText("No file chosen");
		} catch (Exception e) {
			System.err.println("There was a problem opening "+initialFilename);
			e.printStackTrace();
		}
		return new VBox(SPACING_2, new Label(title), output);
	}
	
	
	private static double[][] deepClone(double[][] old) {
		double[][] now = new double[old.length][old[0].length];
		for (int i = 0; i < old.length; i ++)
			for (int j = 0; j < old[i].length; j ++)
				now[i][j] = old[i][j];
		return now;
	}
	
	
	private interface Callback {
		void process(File file) throws IOException;
	}
	

	public static void main(String[] args) {
		launch(args);
	}
	
}
