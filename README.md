# MRSt

Monte Carlo simulations of the time resolved magnetic recoil spectrometer.

## Dependencies

The only dependency is Python with MatPlotLib, which I used for plotting stuff. It needs
to be on your PATH.

Oh, also the Apache Commons math package. I tried so hard not to need JAR files, but I
wasn't strong enough. I guess they're math was just better than mine.

`ConsoleEvaluator.java` is designed to be run from the command line. You can do it like
this:
~~~~bash
$ javac -cp [lib]/commons-math3-3.6.1/commons-math3-3.6.1.jar src/*.java
$ java -cp [lib]/commons-math3-3.6.1/commons-math3-3.6.1.jar:src/ app.ConsoleEvaluator [low|med|high] [N]
~~~~
where `lib` is the path from the root directory to that stupid JAR file, the initial
argument is the configuration to evaluate, and `N` is the number of evaluations to do.
