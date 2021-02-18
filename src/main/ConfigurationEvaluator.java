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
import javafx.scene.control.TextField;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import main.CSV.COSYMapping;
import main.MRSt.ErrorMode;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class ConfigurationEvaluator extends Application {
	
	private static final Particle ION = Particle.D;
	private static final File STOPPING_POWER_FILE = new File("data/stopping_power_deuterons.csv");
	private static final double COSY_MINIMUM_ENERGY = 10.7e6;
	private static final double COSY_MAXIMUM_ENERGY = 14.2e6;
	private static final double COSY_REFERENCE_ENERGY = 12.45e6;
	
	private static final int SPACING_0 = 16;
	private static final int SPACING_1 = 10;
	private static final int SPACING_2 = 4;
	
	private static final int NUM_YIELDS = 18;
	
	private static final String[] HEADERS = {
		"Yield factor", "Temperature factor", "Down-scatter factor", "Velocity shift (μm/ns)",
		"Computation time (s)", "Total yield (10^15)", "Bang time (ns)",
		"Burn width (ns)", "Burn skewness", "Burn kurtosis",
		"Max \u03C1R (ns)", "\u03C1R at max (g/cm^2)",
		"Burn-average \u03C1R (g/cm^2)", "d\u03C1R/dt at BT (g/cm^2/ns)",
		"Burn-average Ti (keV)", "dTi/dt at BT (keV/ns)",
		"Burn-average vi (km/s)", "dvi/dt at BT (km/s/ns)",
		"Peak Ti (keV)", "Energy confinement time (ns)",
		"Cooling Ti (keV)", "Cooling time (ns)"};
	private static final String[] HEADERS_WITH_ERRORS = new String[(HEADERS.length-4)*2+4];
	
	private Spinner<Double> foilDistance;
	private Spinner<Double> foilRadius;
	private Spinner<Double> foilThickness;
	private Spinner<Double> apertureDistance;
	private Spinner<Double> apertureWidth;
	private Spinner<Double> apertureHeight;
	private Spinner<Double> focalPlaneTilt;
	private ChoiceBox<Integer> order;
	private CheckBox[] variations;
	private CheckBox errorBars;
	private TextField saveFile;
	
	private double[][] stoppingPowerData;
	private double[][] cosyCoefficients;
	private int[][] cosyExponents;
	private Logger logger;
	
	
	/**
	 * build the GUI and display it.
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public void start(Stage stage) throws NumberFormatException, IOException {
		for (int i = 0; i < HEADERS.length; i ++) {
			if (i < 4)
				HEADERS_WITH_ERRORS[i] = HEADERS[i];
			else {
				HEADERS_WITH_ERRORS[2*(i-4)+4] = HEADERS[i];
				int locationOfTheWordQuoteErrorUnquote = HEADERS[i].indexOf('(') - 1;
				if (locationOfTheWordQuoteErrorUnquote == -2)
					locationOfTheWordQuoteErrorUnquote = HEADERS[i].length();
				HEADERS_WITH_ERRORS[2*(i-4)+5] = HEADERS[i] + " error";
			}
		}
		
		GridPane leftPane = new GridPane();
		leftPane.setHgap(SPACING_1);
		leftPane.setVgap(SPACING_1);
		int row = 0;
		
		this.foilDistance = new Spinner<Double>(0.5, 10.0, 3.0, 0.1);
		foilDistance.setEditable(true);
		leftPane.add(new Label("Foil distance"), 0, row);
		leftPane.add(foilDistance, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilRadius = new Spinner<Double>(0.1, 1.0, 0.4, 0.01);
		foilRadius.setEditable(true);
		leftPane.add(new Label("Foil radius"), 0, row);
		leftPane.add(foilRadius, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilThickness = new Spinner<Double>(5., 500., 100., 5.);
		foilThickness.setEditable(true);
		leftPane.add(new Label("Foil thickness"), 0, row);
		leftPane.add(foilThickness, 1, row);
		leftPane.add(new Label("μm"), 2, row);
		row ++;
		
		this.apertureDistance = new Spinner<Double>(1.0, 10.0, 6.0, 0.5);
		apertureDistance.setEditable(true);
		leftPane.add(new Label("Aper. distance"), 0, row);
		leftPane.add(apertureDistance, 1, row);
		leftPane.add(new Label("m"), 2, row);
		row ++;
		
		this.apertureWidth = new Spinner<Double>(1.0, 50.0, 5.0, 1.0);
		apertureWidth.setEditable(true);
		leftPane.add(new Label("Aper. width"), 0, row);
		leftPane.add(apertureWidth, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.apertureHeight = new Spinner<Double>(1.0, 50.0, 20.0, 1.0);
		apertureHeight.setEditable(true);
		leftPane.add(new Label("Aper. height"), 0, row);
		leftPane.add(apertureHeight, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.focalPlaneTilt = new Spinner<Double>(0.0, 89.9, 66.586, 5.0);
		focalPlaneTilt.setEditable(true);
		leftPane.add(new Label("F. plane angle"), 0, row);
		leftPane.add(focalPlaneTilt, 1, row);
		leftPane.add(new Label("°"), 2, row);
		row ++;
		
		this.order = new ChoiceBox<Integer>(FXCollections.observableArrayList(1, 2, 3));
		order.setValue(3);
		leftPane.add(new Label("Order"), 0, row);
		leftPane.add(order, 1, row);
		row ++;
		
		leftPane.add(chooseFileWidget("COSY map file:", stage, "MRSt_IRF_FP tilted_final.txt",
				(file) -> {
					COSYMapping map = CSV.readCosyCoefficients(file, order.getValue());
					this.cosyCoefficients = map.coefficients;
					this.cosyExponents = map.exponents;
				}), 0, row, 3, 1);
		row ++;
		
		VBox rightPane = new VBox(SPACING_1);
		
		this.stoppingPowerData = CSV.read(STOPPING_POWER_FILE, ',');
		
		this.variations = new CheckBox[4];
		this.variations[0] = new CheckBox("Vary yield");
		this.variations[1] = new CheckBox("Vary temperature");
		this.variations[2] = new CheckBox("Vary density");
		this.variations[3] = new CheckBox("Vary velocity");
		for (CheckBox checkBox: variations) {
			checkBox.setSelected(false);
			rightPane.getChildren().add(checkBox);
		}
		this.variations[0].setSelected(true);
		
		this.errorBars = new CheckBox("Compute error bars");
		this.errorBars.setSelected(true);
		rightPane.getChildren().add(errorBars);
		
		this.saveFile = new TextField("ensemble.csv");
		rightPane.getChildren().add(new VBox(SPACING_2,
				new Label("Output file:"),
				saveFile));
		
		Button execute = new Button("Evaluate!");
		execute.setOnAction((event) -> {
			if (cosyCoefficients == null)
				logger.severe("Please select a COSY map file.");
			else {
				new Thread(() -> {
					MRSt mc = null;
					try {
						mc = new MRSt(
								ION,
								foilDistance.getValue()*1e-3,
								foilRadius.getValue()*1e-3,
								foilThickness.getValue()*1e-6,
								stoppingPowerData,
								apertureDistance.getValue()*1e0,
								apertureWidth.getValue()*1e-3,
								apertureHeight.getValue()*1e-3,
								COSY_MINIMUM_ENERGY,
								COSY_MAXIMUM_ENERGY,
								COSY_REFERENCE_ENERGY,
								cosyCoefficients,
								cosyExponents,
								focalPlaneTilt.getValue(),
								1,
								logger); // make the simulation
					} catch (Exception e) {
						logger.log(Level.SEVERE, e.getMessage(), e);
					}
					
					double[][] results = new double[NUM_YIELDS][HEADERS_WITH_ERRORS.length];
					for (int k = 0; k < NUM_YIELDS; k ++) {
						double[] eBins = null, tBins = null;
						double[][] spec = null;
						try {
							eBins = CSV.readColumn(new File("data/Energy bins.txt"));
							tBins = CSV.readColumn(new File("data/Time bins.txt"));
							spec = CSV.read(new File("data/spectrum.txt"), '\t');
							if (spec.length != eBins.length-1 || spec[0].length != tBins.length-1) {
								logger.info("interpreting a weird spectrum file...");
								spec = MRSt.interpretSpectrumFile(tBins, eBins, spec);
							}
						} catch (ArrayIndexOutOfBoundsException | NumberFormatException e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
						} catch (IOException e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
						}
						
						double yield = (variations[0].isSelected()) ? Math.pow(10, -3.*Math.random()) : 1;
						double temp =  (variations[1].isSelected()) ? Math.exp(2*Math.random() - 1) : 1; // roll the dies on the spectrum modifications
						double downS = (variations[2].isSelected()) ? Math.exp(2*Math.random() - 1) : 1;
						double flow =  (variations[3].isSelected()) ? 200*Math.random()*(2*Math.random() - 1) : 0;
						MRSt.modifySpectrum(tBins, eBins, spec, yield, temp, downS, flow);
						
						ErrorMode errorBars = this.errorBars.isSelected() ?
								ErrorMode.HESSIAN :
								ErrorMode.STATISTICS;
						
						logger.info(String.format("Yn = %f (%d/%d)", yield, k, NUM_YIELDS));
						
						double[] result;
						try {
							result = mc.respond(
									eBins,
									tBins,
									spec,
									errorBars); // and run it many times!
						} catch (Exception e) {
							logger.log(Level.SEVERE, e.getMessage(), e);
							result = null;
						}
						results[k][0] = yield;
						results[k][1] = temp;
						results[k][2] = downS;
						results[k][3] = flow;
						if (result != null)
							System.arraycopy(result, 0, results[k], 4, result.length);
						else
							for (int i = 4; i < results[k].length; i ++)
								results[k][i] = Double.NaN;
						
						if (k%6 == 5 || k == NUM_YIELDS - 1) {
							try {
								CSV.write(results, new File("working/"+saveFile.getText()), ',',
										HEADERS_WITH_ERRORS);
							} catch (IOException e) {
								logger.log(Level.SEVERE, e.getMessage(), e);
							}
							logger.info("Saved ensemble results to working/"+saveFile.getText());
						}
					}
					try {
						ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/view_ensemble.py", saveFile.getText());
						plotPB.start();
					} catch (IOException e) {
						System.err.println("Could not access plotting files and/or scripts.");
					}
				}).start();
			}
		});
		rightPane.getChildren().add(execute);
		
		final TextArea console = new TextArea();
		console.setEditable(false);
		console.setPrefWidth(400);
		console.setFont(Font.font("Monospace"));
		logger = Logger.getLogger("main");
		logger.setLevel(Level.ALL);
		StreamHandler consoleHandler = new StreamHandler() {
			public void publish(LogRecord record) {
				Platform.runLater(() -> {
					console.appendText(String.format("%7s: %s\n",
							record.getLevel().toString(), record.getMessage()));
				});
			}
		};
		consoleHandler.setLevel(Level.FINER);
		logger.addHandler(consoleHandler);
		rightPane.getChildren().add(console);
		
		StackPane root = new StackPane();
		root.setPadding(new Insets(SPACING_0));
		root.getChildren().add(new HBox(SPACING_0, leftPane, rightPane));
		Scene scene = new Scene(root);
		
		stage.setTitle("Configuration evaluator");
		stage.setScene(scene);
		stage.show();
	}
	
	
	/**
	 * create a consistent-looking file selection thingy with a button and label, and bind it
	 * to the given File.
	 * @param title the title of the window
	 * @param stage
	 * @param initialFilename the filename of the file that it will try to load initially
	 * @param action the action to take with the file once it's chosen
	 * @return
	 */
	private static Region chooseFileWidget(String title, Stage stage, String initialFilename,
			Callback action) {
		Label label = new Label();
		Button button = new Button("Chose file…");
		FileChooser fileChooser = new FileChooser();
		fileChooser.setTitle("Load "+title);
		fileChooser.setInitialDirectory(new File("data/"));
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
			action.process(new File("data/"+initialFilename));
		} catch (IOException e) {
			label.setText("No file chosen");
		} catch (Exception e) {
			System.err.println("There was a problem opening "+initialFilename);
		}
		return new VBox(SPACING_2, new Label(title), output);
	}
	
	
	private static interface Callback {
		void process(File file) throws IOException;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		launch(args);
	}
	
}
