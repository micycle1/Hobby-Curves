package micycle.hobbycurves;

import org.apache.commons.math3.complex.Complex;

/**
 * Models a point on the curve, having auxiliary parameters about the curve
 * passing through it.
 */
class Knot {

	/** Knot coordinates. */
	final Complex cmplx;
	/** Tension of curve at this knot. */
	double alpha = 1.0;
	double beta = 1.0;
	/** Distance between this and next knots. */
	double distance = 0.0;
	/** Curve <b>departure</b> angle (radians) at this knot. */
	double theta = 0.0;
	/** Curve <b>arrival</b> angle (radians) at this knot. */
	double phi = 0.0;
	/** Curve <b>turning</b> angle (radians) at this knot. */
	double psi = 0.0;

	Knot(Complex complex, double alpha, double beta) {
		this.cmplx = complex;
		this.alpha = alpha;
		this.beta = beta;
	}

}