package micycle.hobbycurves;

import org.apache.commons.math3.complex.Complex;

class HobbyPoint {
	
	public HobbyPoint(Complex complex, double alpha, double beta) {
		this.cmplx = complex;
		this.alpha = alpha;
		this.beta = beta;
	}

	// Algorithm parameters.
	Complex cmplx = new Complex(0, 0);
// In what follows, we use Knuth's notation in our variable names
	double alpha = 1.0;
	double beta = 1.0;
	double d_val = 0.0; // Distance between this point and next.
	double theta = 0.0; // Angle of polygonal line from this point to next.
	double phi = 0.0; // Offset angle.
	double psi = 0.0; // Another offset angle
}