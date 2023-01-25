package micycle.hobbycurves;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * Class for generating smooth interpolating splines using the algorithm
 * described in the paper "Smooth, Easy to Compute Interpolating Splines" by
 * John D. Hobby.
 * <p>
 * Hobby Curves consist of a sequence of quadratic bezier curves which smootly
 * pass through a given sequence of knots.
 * 
 * @author Michael Carleton
 * @author Luke Trujillo
 *
 */
public class HobbyCurve {

	// implements https://github.com/ltrujello/Hobby_Curve_Algorithm/

	private List<Knot> knots;
	private double[][] inputPts;
	private double[][] ctrlPts;
	private final boolean closed;
	private final double beginCurl;
	private final double endCurl;
	private final int numPoints;

	/**
	 * Constructor for creating a new "Hobby Curve".
	 *
	 * @param inputPoints a 2D array of [x, y] coordinates of the data knots that
	 *                    the curve should pass through
	 * @param tension     a parameter that controls the tension of the curve's
	 *                    "knots". A value of 1 is a good starting point.
	 * @param closed      a boolean value indicating whether the curve should form a
	 *                    closed loop, joining the first and last knots
	 * @param beginCurl   a value that controls the amount of curl at the start of
	 *                    the curve (the bezier originating at the starting point).
	 *                    Accepts any value, but should generally be positive; a
	 *                    value of 1 is a good starting point. This value is only
	 *                    effective if the curve is unclosed.
	 * @param endCurl     a value that controls the amount of curl at the end of the
	 *                    curve (the bezier terminating at the final point). Accepts
	 *                    any value, but should generally be positive; a value of 1
	 *                    is a good starting point. This value is only effective if
	 *                    the curve is unclosed.
	 */
	public HobbyCurve(double[][] inputPoints, double tension, boolean closed, double beginCurl, double endCurl) {
		if (inputPoints.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two knots");
		}
		this.knots = new ArrayList<>(inputPoints.length);
		this.inputPts = inputPoints;
		for (double[] point : inputPoints) {
			Knot new_point = new Knot(new Complex(point[0], point[1]), 1.0 / tension, 1.0 / tension);
			knots.add(new_point);
		}
		this.closed = closed;
		this.beginCurl = beginCurl;
		this.endCurl = endCurl;
		this.numPoints = knots.size();
	}

	/**
	 * Returns the 4 parameters (each an (x,y) coordinate pair) of every quadratic
	 * bezier curve comprising this Hobby Curve. Each row of the output has the
	 * form...
	 * <p>
	 * <code>[ap1.x, ap1.y, cp1.x ,cp1.y, cp2.x, cp2.y, ap2.x, ap2.y]</code>
	 * <p>
	 * ...where <code>ap</code> denotes "anchor point" and <code>cp</code> denotes
	 * "control point".
	 * 
	 * @return a list of bezier curve parameters
	 */
	public double[][] getBeziers() {
		getBezierCtrlPts();
		double[][] bezierParams = new double[numPoints - (closed ? 0 : 1)][8];
		for (int i = 0; i < bezierParams.length; i++) { // bezier between ith and i+1th point
			double[] ap1 = inputPts[i]; // anchor point 1
			double[] cp1 = ctrlPts[2 * i]; // control point 1
			double[] cp2 = ctrlPts[2 * i + 1];
			double[] ap2 = inputPts[(i + 1) % inputPts.length];
			bezierParams[i][0] = ap1[0]; // x
			bezierParams[i][1] = ap1[1]; // y
			bezierParams[i][2] = cp1[0];
			bezierParams[i][3] = cp1[1];
			bezierParams[i][4] = cp2[0];
			bezierParams[i][5] = cp2[1];
			bezierParams[i][6] = ap2[0];
			bezierParams[i][7] = ap2[1];
		}
		return bezierParams;
	}

	/**
	 * Calculates the control points for each of the quadratic bezier curves
	 * comprising this Hobby Curve. In the output, each sequential pair of (x, y)
	 * coordinates denote the two control points of a single bezier.
	 * 
	 * @return 2d array of control point (x, y) coordinates
	 */
	public double[][] getBezierCtrlPts() {
		if (ctrlPts == null) {
			calculateDVals();
			calculatePsiVals();
			calculateThetaVals();
			calculatePhiVals();
			calculateCtrlPts();
		}
		return ctrlPts;
	}

	/** Calculates the pairwise distances between the knots. */
	private void calculateDVals() {
		// Skip last point if path is non-closed
		int end = closed ? numPoints : numPoints - 1;
		for (int i = 0; i < end; i++) {
			Knot z_i = knots.get(i);
			Knot z_j = knots.get((i + 1) % numPoints);
			z_i.distance = z_i.cmplx.subtract(z_j.cmplx).abs();
		}
	}

	/** Calculates the psi values by subtracting pairwise phases. */
	private void calculatePsiVals() {
		// Skip first and last point if path is non-closed
		int start = closed ? 0 : 1;
		int end = closed ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			Knot z_h = (i == 0 ? knots.get(numPoints - 1) : knots.get(i - 1));
			Knot z_i = knots.get(i);
			Knot z_j = knots.get((i + 1) % numPoints);
			Complex polygonal_turn = z_j.cmplx.subtract(z_i.cmplx).divide(z_i.cmplx.subtract(z_h.cmplx));
			z_i.psi = Math.atan2(polygonal_turn.getImaginary(), polygonal_turn.getReal());
		}
	}

	/**
	 * Calculates the theta values by creating a linear system whose solutions are
	 * the values.
	 */
	private void calculateThetaVals() {
		final RealVector A = new ArrayRealVector(numPoints);
		final RealVector B = new ArrayRealVector(numPoints);
		final RealVector C = new ArrayRealVector(numPoints);
		final RealVector D = new ArrayRealVector(numPoints);
		final RealVector R = new ArrayRealVector(numPoints);

		// Calculate the entries of the five vectors.
		// Skip first and last point if path is non-closed.
		int start = closed ? 0 : 1;
		int end = closed ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			Knot z_h = (i == 0 ? knots.get(numPoints - 1) : knots.get(i - 1));
			Knot z_i = knots.get(i);
			Knot z_j = knots.get((i + 1) % numPoints);

			A.setEntry(i, z_h.alpha / (z_i.beta * z_i.beta * z_h.distance));
			B.setEntry(i, (3 - z_h.alpha) / (z_i.beta * z_i.beta * z_h.distance));
			C.setEntry(i, (3 - z_j.beta) / (z_i.alpha * z_i.alpha * z_i.distance));
			D.setEntry(i, z_j.beta / (z_i.alpha * z_i.alpha * z_i.distance));
			R.setEntry(i, -B.getEntry(i) * z_i.psi - D.getEntry(i) * z_j.psi);
		}

		RealMatrix M = new Array2DRowRealMatrix(numPoints, numPoints);
		// Set up matrix M such that the soln. Mx = R are the theta values.
		for (int i = start; i < end; i++) {
			// Fill i-th row of M
			if (i == 0) {
				M.setEntry(i, numPoints - 1, A.getEntry(i));
			} else {
				M.setEntry(i, i - 1, A.getEntry(i));
			}
			M.setEntry(i, i, B.getEntry(i) + C.getEntry(i));
			M.setEntry(i, (i + 1) % numPoints, D.getEntry(i));
		}

		// Special formulas for first and last rows of M with non-closed paths.
		if (!closed) {
			// First row of M
			double alpha_0 = knots.get(0).alpha;
			double beta_1 = knots.get(1).beta;
			double xi_0 = (alpha_0 * alpha_0 * beginCurl) / (beta_1 * beta_1);
			M.setEntry(0, 0, alpha_0 * xi_0 + 3 - beta_1);
			M.setEntry(0, 1, (3 - alpha_0) * xi_0 + beta_1);
			R.setEntry(0, -((3 - alpha_0) * xi_0 + beta_1) * knots.get(1).psi);
			// Last row of M
			double alpha_n_1 = knots.get(numPoints - 2).alpha;
			double beta_n = knots.get(numPoints - 1).beta;
			double xi_n = (beta_n * beta_n * endCurl) / (alpha_n_1 * alpha_n_1);
			M.setEntry(numPoints - 1, numPoints - 2, (3 - beta_n) * xi_n + alpha_n_1);
			M.setEntry(numPoints - 1, numPoints - 1, (beta_n * xi_n + 3 - alpha_n_1));
			R.setEntry(numPoints - 1, 0);
		}

		// Solve for theta values
		DecompositionSolver solver = new LUDecomposition(M).getSolver();
		RealVector thetas = solver.solve(R);
		// Assign theta values to each Knot
		for (int i = 0; i < numPoints; i++) {
			knots.get(i).theta = thetas.getEntry(i);
		}

	}

	/**
	 * Calculates the phi_k values via the relationship
	 * <code>theta_k + phi_k + psi_k = 0</code>.
	 */
	private void calculatePhiVals() {
		for (Knot point : knots) {
			point.phi = -(point.psi + point.theta);
		}
	}

	/** Calculates the Bezier control points from z_i to z_{i+1}. */
	private void calculateCtrlPts() {
		int end = closed ? numPoints : numPoints - 1; // Skip last point if path is non-closed
		ctrlPts = new double[end * 2][2];
		for (int i = 0; i < end; i++) {
			final Knot z_i = knots.get(i);
			final Knot z_j = knots.get((i + 1) % numPoints);
			final double rho_coefficient = z_i.alpha * velocity(z_i.theta, z_j.phi);
			final double sigma_coefficient = z_j.beta * velocity(z_j.phi, z_i.theta);

			Complex ctrl_point_a = z_i.cmplx.add(z_j.cmplx.subtract(z_i.cmplx).multiply(rho_coefficient / 3)
					.multiply(new Complex(Math.cos(z_i.theta), Math.sin(z_i.theta))));
			Complex ctrl_point_b = z_j.cmplx.subtract(z_j.cmplx.subtract(z_i.cmplx).multiply(sigma_coefficient / 3)
					.multiply(new Complex(Math.cos(-z_j.phi), Math.sin(-z_j.phi))));

			ctrlPts[i * 2] = new double[] { ctrl_point_a.getReal(), ctrl_point_a.getImaginary() };
			ctrlPts[i * 2 + 1] = new double[] { ctrl_point_b.getReal(), ctrl_point_b.getImaginary() };
		}
	}

	/** Metafont's velocity function. */
	private static double velocity(final double theta, final double phi) {
		final double sinTheta = Math.sin(theta);
		final double sinPhi = Math.sin(phi);
		final double cosTheta = Math.cos(theta);
		final double cosPhi = Math.cos(phi);
		double numerator = 2 + Math.sqrt(2) * (sinTheta - (1. / 16) * sinPhi) * (sinPhi - (1. / 16) * sinTheta) * (cosTheta - cosPhi);
		double denominator = (1 + (1. / 2) * (Math.sqrt(5) - 1) * cosTheta + (1. / 2) * (3 - Math.sqrt(5)) * cosPhi);
		return numerator / denominator;
	}

}
