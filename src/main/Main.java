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
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.ScrollPane;
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
	
	public static final int WIDTH = 800;
	public static final int HEIGHT = 600;
	
	public static final int N_ROWS = 56;
	public static final int M_COLS = 5;
	
	
	public void start(Stage stage) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader("data/MRSt_IRF_FP tilted.txt"));
		double[][] matrix = new double[N_ROWS][M_COLS];
		for (int i = 0; i < N_ROWS; i ++) {
			String[] parts = in.readLine().split("\\s+");
			for (int j = 0; j < M_COLS; j ++)
				matrix[i][j] = Double.parseDouble(parts[j+1]);
		}
		in.close();
		
		Text text = new Text(WIDTH, HEIGHT, Arrays.deepToString(matrix).replace("], [", "]\n["));
		
		ScrollPane scroll = new ScrollPane(text);
		
		StackPane root = new StackPane();
		root.getChildren().add(new VBox(5, scroll, new Button("button")));
		Scene scene = new Scene(root, WIDTH, HEIGHT);
		
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
