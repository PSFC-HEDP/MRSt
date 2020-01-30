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
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class Main extends Application {
	
	private static final Particle ION = Particle.D;
	private static final File STOPPING_POWER_FILE = new File("data/stopping_power_deuterons.csv");
	private static final double COSY_MINIMUM_ENERGY = 10.7e6;
	private static final double COSY_MAXIMUM_ENERGY = 14.2e6;
	private static final double COSY_REFERENCE_ENERGY = 12.45e6;
	private static final int NUM_BINS = 150;
	private static final int SPACING_0 = 16;
	private static final int SPACING_1 = 10;
	private static final int SPACING_2 = 4;
	
	private Spinner<Double> foilDistance;
	private Spinner<Double> foilRadius;
	private Spinner<Double> foilThickness;
	private Spinner<Double> apertureDistance;
	private Spinner<Double> apertureWidth;
	private Spinner<Double> apertureHeight;
	private Spinner<Double> focalPlaneTilt;
	private ChoiceBox<Integer> order;
	private double[][] stoppingPowerData;
	private double[][] cosyCoefficients;
	private int[][] cosyExponents;
	private double[] timeBins;
	private double[] energyBins;
	private double[][] spectrum;
	private Logger logger;
	
	
	/**
	 * build the GUI and display it.
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public void start(Stage stage) throws NumberFormatException, IOException {
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
		
		this.foilRadius = new Spinner<Double>(0.1, 1.0, 0.3, 0.01); // TODO maybe throw a warning if the radius >~ the distance
		foilRadius.setEditable(true);
		leftPane.add(new Label("Foil radius"), 0, row);
		leftPane.add(foilRadius, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilThickness = new Spinner<Double>(10., 500., 80., 5.);
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
		
		this.apertureWidth = new Spinner<Double>(1.0, 50.0, 4.0, 1.0);
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
		
		this.focalPlaneTilt = new Spinner<Double>(0.0, 89.9, 70.3, 5.0);
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
		
		leftPane.add(chooseFileWidget("COSY data file:", stage, (file) -> {
			this.cosyCoefficients = CSV.readCosyCoefficients(file, order.getValue());
			this.cosyExponents = CSV.readCosyExponents(file, order.getValue());
		}), 0, row, 3, 1);
		row ++;
		
		VBox rightPane = new VBox(SPACING_1);
		
		rightPane.getChildren().add(chooseFileWidget("Energy bin file:", stage, (file) -> {
			this.energyBins = CSV.readColumn(file);
		}));
		row ++;
		
		rightPane.getChildren().add(chooseFileWidget("Time bin file:", stage, (file) -> {
			this.timeBins = CSV.readColumn(file);
		}));
		row ++;
		
		rightPane.getChildren().add(chooseFileWidget("Spectrum file:", stage, (file) -> {
			this.spectrum = CSV.read(file, '\t');
			this.spectrum = MonteCarlo.correctSpectrum(timeBins, energyBins, spectrum); // TODO: move this part to a separate callback where it is also plotted
		}));
		row ++;
		
		this.stoppingPowerData = CSV.read(STOPPING_POWER_FILE, ',');
		
		Button execute = new Button("Compute!");
		execute.setOnAction((event) -> {
			new Thread(() -> {
				MonteCarlo mc = new MonteCarlo(
						ION, foilDistance.getValue()*1e-3, foilRadius.getValue()*1e-3,
						foilThickness.getValue()*1e-6, stoppingPowerData, apertureDistance.getValue()*1e0,
						apertureWidth.getValue()*1e-3, apertureHeight.getValue()*1e-3,
						COSY_MINIMUM_ENERGY, COSY_MAXIMUM_ENERGY, COSY_REFERENCE_ENERGY,
						cosyCoefficients, cosyExponents, focalPlaneTilt.getValue(), NUM_BINS, logger);
				mc.respond(energyBins, timeBins, spectrum);
				try {
					plotHeatmap(timeBins, energyBins, spectrum,
							"Time (ns)", "Energy (MeV)", "Spectrum", logger);
					plotHeatmap(mc.getMeasuredTimeBins(), mc.getPositionBins(), mc.getMeasuredSpectrum(),
							"Time (ns)", "Position (cm)", "Response", logger);
					plotHeatmap(mc.getInferredTimeBins(), mc.getEnergyBins(), mc.getInferredSpectrum(),
							"Time (ns)", "Energy (MeV)", "Inferred", logger);
				} catch (IOException e) {
					logger.severe(e.getMessage());
				}
			}).start();
		});
		rightPane.getChildren().add(execute);
		
		final TextArea console = new TextArea();
		console.setEditable(false);
		console.setPrefWidth(400);
		logger = Logger.getLogger("main");
		logger.addHandler(new StreamHandler() {
			public void publish(LogRecord record) {
				Platform.runLater(() -> {
					console.appendText(String.format("%7s: %s\n",
							record.getLevel().toString(), record.getMessage()));
				});
			}
		});
		rightPane.getChildren().add(console);
		
		StackPane root = new StackPane();
		root.setPadding(new Insets(SPACING_0));
		root.getChildren().add(new HBox(SPACING_0, leftPane, rightPane));
		Scene scene = new Scene(root);
		
		stage.setTitle("MRSt");
		stage.setScene(scene);
		stage.show();
	}
	
	
	private static void plotHeatmap(double[] x, double[] y, double[][] z,
			String xlabel, String ylabel, String title, Logger logger) throws IOException {
		CSV.writeColumn(x, new File(String.format("working/%s_x.csv", title)));
		CSV.writeColumn(y, new File(String.format("working/%s_y.csv", title)));
		CSV.write(z, new File(String.format("working/%s_z.csv", title)), ',');
		ProcessBuilder plotPB = new ProcessBuilder("python", "src/python/plot.py", xlabel, ylabel, title);
		plotPB.start();
	}
	
	
	/**
	 * create a consistent-looking file selection thingy with a button and label, and bind it
	 * to the given File.
	 * @param name
	 * @param bound
	 * @return
	 */
	private static Region chooseFileWidget(String title, Stage stage, Callback action) {
		Label label = new Label("No file chosen.");
		Button button = new Button("Chose file…");
		FileChooser fileChooser = new FileChooser();
		fileChooser.setTitle("Load data file");
		fileChooser.setInitialDirectory(new File("data/"));
		fileChooser.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("Data files", "*.csv", "*.txt"),
				new FileChooser.ExtensionFilter("All files", "*.*"));
		button.setOnAction((event) -> {
			final File chosen = fileChooser.showOpenDialog(stage);
			if (chosen != null) {
				label.setText(chosen.getName());
				try {
					action.process(chosen);
				} catch (IOException e) {
					System.err.println("Could not open file."); // TODO remind myself how Alerts work
				} catch (Exception e) {
					System.err.println("There was a problem opening the file.");
					e.printStackTrace();
				}
			}
		});
		HBox output = new HBox(SPACING_2, button, label);
		label.setMaxWidth(250 - button.getPrefWidth());
		output.setPrefWidth(250);
		output.setAlignment(Pos.CENTER_LEFT);
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
