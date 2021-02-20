/**
 * 
 */
package main;

/**
 * some equations about implosion dynamics
 * @author Justin Kunimune
 */
public class ImplosionModelling {

//	private static final double eV = Math.abs(Particle.E.charge);
//	private static final double keV = 1e3*Math.abs(Particle.E.charge);
//	private static final double MeV = 1e6*Math.abs(Particle.E.charge);
//	private static final double μm = 1e-6;
//	private static final double ns = 1e-9;
//	private static final double ps = 1e-12;
//	private static final double GPa = 1e9;
//	
//	
//	
//	/**
//	 * solve that transcendental equation at each time step to get the pressure as a function
//	 * of time, given the initial hot spot conditions and some other stuff.
//	 * @param yield neutron yield [10^15/ns]
//	 * @param temperature ion temperature [keV]
//	 * @param confinementTime energy confinement time [ns]
//	 * @param peakTemperature maximum ion temperature [keV]
//	 * @param initPressure initial hot spot pressure [GPa]
//	 * @param initVolume initial hot spot volume [m^3]
//	 * @return
//	 */
//	public static Quantity[] solveForPressure(Quantity[] yield, Quantity[] temperature,
//			Quantity confinementTime, Quantity peakTemperature, Quantity initPressure, Quantity initVolume) {
//		Quantity[] pressure = new Quantity[yield.length];
//		double τECT = confinementTime.times(ns).value, TPeak = peakTemperature.times(keV).value;
//		double p0 = initPressure.times(GPa).value, V0 = initVolume.value;
//		double H = H(peakTemperature.value);
//		for (int i = 0; i < pressure.length; i ++) {
//			double Yn = yield[i].times(1e15/ns).value, Ti = temperature[i].times(keV).value;
//			double σv = reactivity(temperature[i]).value;
//			double leftHandSide = 16*Yn*Ti/σv;
//			double min = 0, max = 2*Math.pow(leftHandSide/V0/Math.pow(p0, 3/5.), 5/7.);
//		}
//		return pressure;
//	}
//	
//	/**
//	 * compute the volume as a function of yield, temperature, and pressure, using the ideal
//	 * gas law and DT reactivity function
//	 * @param yield neutron yield [10^15/ns]
//	 * @param temperature ion temperature [keV]
//	 * @param pressure hot spot pressure [GPa]
//	 * @return hot spot pressure [m^3]
//	 */
//	public static Quantity[] solveForVolume(Quantity[] yield, Quantity[] temperature, Quantity[] pressure) {
//		Quantity[] volume = new Quantity[yield.length];
//		for (int i = 0; i < volume.length; i ++) {
//			Quantity Yn = yield[i].times(1e15/ns), Ti = temperature[i].times(keV), p = pressure[i].times(GPa);
//			Quantity σv = reactivity(temperature[i]);
//			volume[i] = Yn.over(p.pow(2)).times(Ti).over(σv);
//		}
//		return volume;
//	}
//	
//	public static Quantity solveForIHF(Quantity[] temperature, Quantity[] volume) {
//		
//	}
//	
//	public static double H(double T) {
//		
//	}
//	
//	/**
//	 * look up the ⟨σv⟩ reactivity coefficient for temperature
//	 * @param temperature ion temperature [keV]
//	 * @return mean cross section velocity product [m^3/s]
//	 */
//	public static Quantity reactivity(Quantity temperature) {
//		return temperature.pow(-1/3.).times(-19.94).exp().times(temperature.pow(-2/3.)).times(3.68e-18);
//	}

}
