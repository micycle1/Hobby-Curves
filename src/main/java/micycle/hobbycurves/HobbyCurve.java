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
 * 
 * @author Michael Carleton
 *
 */
public class HobbyCurve {

	// implements https://github.com/ltrujello/Hobby_Curve_Algorithm/

	private List<HobbyPoint> points;
	private double[][] ctrlPts;
	private final boolean cyclic;
	private final double beginCurl;
	private final double endCurl;
	private final int numPoints;

	/**
	 * Constructor for creating a new interpolating spline ("Hobby Curve").
	 *
	 * @param inputPoints a 2D array of x,y coordinates of the data points that the
	 *                    curve should approximate
	 * @param tension     a value that controls the tension of the curve "knots", it
	 *                    should be between 0 and 1, with 0 giving a linear
	 *                    interpolation and 1 giving a tight fit to the data points.
	 *                    // TODO check
	 * @param cyclic      a boolean value indicating whether the spline should be
	 *                    cyclical / form a closed loop.
	 * @param beginCurl   a value that controls the amount of curl at the start of
	 *                    the curve
	 * @param endCurl     a value that controls the amount of curl at the end of the
	 *                    curve
	 */
	public HobbyCurve(double[][] inputPoints, double tension, boolean cyclic, double beginCurl, double endCurl) {
		if (inputPoints.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two points");
		}
		this.points = new ArrayList<>(inputPoints.length);
		for (double[] point : inputPoints) {
			HobbyPoint new_point = new HobbyPoint(new Complex(point[0], point[1]), 1.0 / tension, 1.0 / tension);
			points.add(new_point);
		}
		this.cyclic = cyclic;
		this.beginCurl = beginCurl;
		this.endCurl = endCurl;
		this.numPoints = points.size();
	}

	/**
	 * Calculates the control points of a sequence of Bezier splines which pass
	 * through the input points.
	 * 
	 * @return
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

	/** Calculates the pairwise distances between the points. */
	private void calculateDVals() {
		// Skip last point if path is non-cyclic
		int end = cyclic ? numPoints : numPoints - 1;
		for (int i = 0; i < end; i++) {
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % numPoints);
			z_i.dVal = z_i.cmplx.subtract(z_j.cmplx).abs();
		}
	}

	/** Calculates the psi values by subtracting pairwise phases. */
	private void calculatePsiVals() {
		// Skip first and last point if path is non-cyclic
		int start = cyclic ? 0 : 1;
		int end = cyclic ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			HobbyPoint z_h = (i == 0 ? points.get(numPoints - 1) : points.get(i - 1));
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % numPoints);
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
		// Skip first and last point if path is non-cyclic.
		int start = cyclic ? 0 : 1;
		int end = cyclic ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			HobbyPoint z_h = (i == 0 ? points.get(numPoints - 1) : points.get(i - 1));
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % numPoints);

			A.setEntry(i, z_h.alpha / (z_i.beta * z_i.beta * z_h.dVal));
			B.setEntry(i, (3 - z_h.alpha) / (z_i.beta * z_i.beta * z_h.dVal));
			C.setEntry(i, (3 - z_j.beta) / (z_i.alpha * z_i.alpha * z_i.dVal));
			D.setEntry(i, z_j.beta / (z_i.alpha * z_i.alpha * z_i.dVal));
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

		// Special formulas for first and last rows of M with non-cyclic paths.
		if (!cyclic) {
			// First row of M
			double alpha_0 = points.get(0).alpha;
			double beta_1 = points.get(1).beta;
			double xi_0 = (alpha_0 * alpha_0 * beginCurl) / (beta_1 * beta_1);
			M.setEntry(0, 0, alpha_0 * xi_0 + 3 - beta_1);
			M.setEntry(0, 1, (3 - alpha_0) * xi_0 + beta_1);
			R.setEntry(0, -((3 - alpha_0) * xi_0 + beta_1) * points.get(1).psi);
			// Last row of M
			double alpha_n_1 = points.get(numPoints - 2).alpha;
			double beta_n = points.get(numPoints - 1).beta;
			double xi_n = (beta_n * beta_n * endCurl) / (alpha_n_1 * alpha_n_1);
			M.setEntry(numPoints - 1, numPoints - 2, (3 - beta_n) * xi_n + alpha_n_1);
			M.setEntry(numPoints - 1, numPoints - 1, (beta_n * xi_n + 3 - alpha_n_1));
			R.setEntry(numPoints - 1, 0);
		}

		// Solve for theta values
		DecompositionSolver solver = new LUDecomposition(M).getSolver();
		RealVector thetas = solver.solve(R);
		// Assign theta values to each HobbyPoint
		for (int i = 0; i < numPoints; i++) {
			points.get(i).theta = thetas.getEntry(i);
		}

	}

	/**
	 * Calculates the phi_k values via the relationship
	 * <code>theta_k + phi_k + psi_k = 0</code>.
	 */
	private void calculatePhiVals() {
		for (HobbyPoint point : points) {
			point.phi = -(point.psi + point.theta);
		}
	}

	/** Calculates the Bezier control points from z_i to z_{i+1}. */
	private void calculateCtrlPts() {
		int end = cyclic ? numPoints : numPoints - 1; // Skip last point if path is non-cyclic
		ctrlPts = new double[(numPoints - 1) * 2][2];
		for (int i = 0; i < end; i++) {
			final HobbyPoint z_i = points.get(i);
			final HobbyPoint z_j = points.get((i + 1) % numPoints);
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
