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
import java.io.FileReader;
import java.util.Arrays;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Spinner;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.Stage;


/**
 * the class that handles the GUI.
 * 
 * @author Justin Kunimune
 */
public class Main extends Application {
	
	public static final int N_ROWS = 56;
	public static final int M_COLS = 5;
	
	private static final Particle ION = Particle.D;
	
	private Spinner<Double> foilDistance;
	private Spinner<Double> foilRadius;
	private Spinner<Double> foilThickness;
	private Spinner<Double> apertureDistance;
	private Spinner<Double> apertureWidth;
	private Spinner<Double> apertureHeight;
	private Label cosyFile;
	private Label timeBinFile;
	private Label energyBinFile;
	private Label spectrumFile;
	
	private double[][] cosyCoefficients;
	private int[][] cosyExponents;
	
	
	public void start(Stage stage) throws Exception {
		GridPane leftPane = new GridPane();
		leftPane.setHgap(6);
		leftPane.setVgap(6);
		int row = 0;
		
		this.foilDistance = new Spinner<Double>(0.5, 10.0, 3.0);
		leftPane.add(new Label("Foil distance"), 0, row);
		leftPane.add(foilDistance, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilRadius = new Spinner<Double>(0.1, 1.0, 0.3); // TODO maybe throw a warning if the radius >~ the distance
		leftPane.add(new Label("Foil radius"), 0, row);
		leftPane.add(foilRadius, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.foilThickness = new Spinner<Double>(10., 500., 80.);
		leftPane.add(new Label("Foil thickness"), 0, row);
		leftPane.add(foilThickness, 1, row);
		leftPane.add(new Label("μm"), 2, row);
		row ++;
		
		this.apertureDistance = new Spinner<Double>(1.0, 10.0, 6.0);
		leftPane.add(new Label("Aper. distance"), 0, row);
		leftPane.add(apertureDistance, 1, row);
		leftPane.add(new Label("m"), 2, row);
		row ++;
		
		this.apertureWidth = new Spinner<Double>(1.0, 50.0, 4.0);
		leftPane.add(new Label("Aper. width"), 0, row);
		leftPane.add(apertureWidth, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		this.apertureHeight = new Spinner<Double>(1.0, 50.0, 20.0);
		leftPane.add(new Label("Aper. height"), 0, row);
		leftPane.add(apertureHeight, 1, row);
		leftPane.add(new Label("mm"), 2, row);
		row ++;
		
		leftPane.add(new Label("COSY matrix file:"), 0, row, 3, 1);
		row ++;
		
		cosyFile = new Label("No file chosen.");
		Button chooseCosyFile = new Button("Choose...");
		leftPane.add(new HBox(3, chooseCosyFile, cosyFile), 0, row, 3, 1);
		row ++;
		
		VBox middlePane = new VBox(6);
		
		GridPane middleSubpane = new GridPane();
		middleSubpane.setHgap(6);
		middleSubpane.setVgap(6);
		middlePane.getChildren().add(middleSubpane);
		row = 0;
		
		middleSubpane.add(new Label("Time bin file:"), 0, row, 3, 1);
		row ++;
		
		timeBinFile = new Label("No file chosen.");
		Button chooseTimeBinFile = new Button("Choose...");
		middleSubpane.add(new HBox(3, chooseTimeBinFile, timeBinFile), 0, row, 3, 1);
		row ++;
		
		middleSubpane.add(new Label("Energy bin file:"), 0, row, 3, 1);
		row ++;
		
		energyBinFile = new Label("No file chosen.");
		Button chooseEnergyBinFile = new Button("Choose...");
		middleSubpane.add(new HBox(3, chooseEnergyBinFile, energyBinFile), 0, row, 3, 1);
		row ++;
		
		middleSubpane.add(new Label("Spectrum file:"), 0, row, 3, 1);
		row ++;
		
		spectrumFile = new Label("No file chosen.");
		Button chooseSpectrumFile = new Button("Choose...");
		middleSubpane.add(new HBox(3, chooseSpectrumFile, spectrumFile), 0, row, 3, 1);
		row ++;
		
		VBox rightPane = new VBox(6);
		
		Button execute = new Button("Compute!");
		rightPane.getChildren().add(execute);
		
		StackPane root = new StackPane();
		root.setPadding(new Insets(12));
		root.getChildren().add(new HBox(12, leftPane, middlePane, rightPane));
		Scene scene = new Scene(root);
		
		stage.setTitle("MRSt");
		stage.setScene(scene);
		stage.show();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		launch(args);
	}
	
}
