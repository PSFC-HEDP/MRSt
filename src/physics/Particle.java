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
package physics;


/**
 * a simple enum to store various particle parameters in a relatively readable way.
 * 
 * @author Justin Kunimune
 */
public enum Particle {
	
	/** electron */
	E(9.10938356e-31, -1.60217662e-19, "electron"),
	/** proton */
	P(1.67262192e-27,  1.60217662e-19, "proton"),
	/** neutron */
	N(1.67492750e-27,  0.0,            "neutron"),
	/** deuteron */
	D(3.34449463e-27,  1.60217662e-19, "deuteron");
	
	public final double mass; // [kg]
	public final double charge; // [C]
	public final String name;
	
	Particle(double mass, double charge, String name) {
		this.mass = mass;
		this.charge = charge;
		this.name = name;
	}
}
