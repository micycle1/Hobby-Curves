package micycle.hobbycurves;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import micycle.hobbycurves.Knot.Coordinate;

/**
 * Class for generating smooth interpolating splines using the algorithm
 * described in the paper "Smooth, Easy to Compute Interpolating Splines" by
 * John D. Hobby.
 * <p>
 * Hobby Curves consist of a sequence of cubic bezier curves which smoothly pass
 * through a given sequence of knots.
 * 
 * @author Michael Carleton
 * @author Luke Trujillo
 *
 */
public class HobbyCurve implements IHobbyCurve {

	// implements https://github.com/ltrujello/Hobby_Curve_Algorithm/
	// https://tex.stackexchange.com/questions/54771/curve-through-a-sequence-of-points-with-metapost-and-tikz

	private List<Knot> knots;
	private double[][] bezierParams;
	private final boolean closed;
	private final double beginCurl;
	private final double endCurl;
	private final int numPoints;

	/**
	 * Constructor for creating a new "Hobby Curve".
	 *
	 * @param inputPoints a 2D array of [x, y] coordinates of the data knots that
	 *                    the curve should pass through.
	 * @param tension     a parameter that controls the tension of the curve's
	 *                    "knots". A value of 1 is a good starting point.
	 * @param closed      a boolean value indicating whether the curve should form a
	 *                    closed loop, joining the first and last knots
	 * @param beginCurl   a value that controls the amount of curl at the start of
	 *                    the curve (the angle of incoming bezier originating at the
	 *                    starting point). Accepts any value, but should generally
	 *                    be positive; a value of 1 is a good starting point. This
	 *                    value is only effective if the curve is unclosed.
	 * @param endCurl     a value that controls the amount of curl at the end of the
	 *                    curve (the angle out outgoing bezier terminating at the
	 *                    final point). Accepts any value, but should generally be
	 *                    positive; a value of 1 is a good starting point. This
	 *                    value is only effective if the curve is unclosed.
	 */
	public HobbyCurve(double[][] inputPoints, double tension, boolean closed, double beginCurl, double endCurl) {
		if (inputPoints.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two knots");
		}
		this.knots = new ArrayList<>(inputPoints.length);
		for (double[] point : inputPoints) {
			Knot knot = new Knot(new Coordinate(point[0], point[1]), 1.0 / tension, 1.0 / tension);
			knots.add(knot);
		}
		this.closed = closed;
		this.beginCurl = beginCurl;
		this.endCurl = endCurl;
		this.numPoints = knots.size();
	}

	/**
	 * Constructor for creating a new "Hobby Curve", where tension is specified per
	 * knot.
	 *
	 * @param inputPoints a 2D array of [x, y] coordinates of the data knots that
	 *                    the curve should pass through.
	 * @param tensions    a list of tension parameters that control the tension of
	 *                    the curve at each corresponding knots. A value of 1 is a
	 *                    good starting point.
	 * @param closed      a boolean value indicating whether the curve should form a
	 *                    closed loop, joining the first and last knots
	 * @param beginCurl   a value that controls the amount of curl at the start of
	 *                    the curve (the angle of incoming bezier originating at the
	 *                    starting point). Accepts any value, but should generally
	 *                    be positive; a value of 1 is a good starting point. This
	 *                    value is only effective if the curve is unclosed.
	 * @param endCurl     a value that controls the amount of curl at the end of the
	 *                    curve (the angle out outgoing bezier terminating at the
	 *                    final point). Accepts any value, but should generally be
	 *                    positive; a value of 1 is a good starting point. This
	 *                    value is only effective if the curve is unclosed.
	 */
	public HobbyCurve(double[][] inputPoints, double[] tensions, boolean closed, double beginCurl, double endCurl) {
		if (inputPoints.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two knots");
		}
		if (inputPoints.length != tensions.length) {
			throw new IllegalArgumentException("Points and tensions arrays are different sizes. They should be equal.");
		}
		this.knots = new ArrayList<>(inputPoints.length);
		for (int i = 0; i < inputPoints.length; i++) {
			double[] point = inputPoints[i];
			double tension = tensions[i];
			Knot knot = new Knot(new Coordinate(point[0], point[1]), 1.0 / tension, 1.0 / tension);
			knots.add(knot);
		}
		this.closed = closed;
		this.beginCurl = beginCurl;
		this.endCurl = endCurl;
		this.numPoints = knots.size();
	}

	/**
	 * Returns the 4 parameters (each an (x,y) coordinate pair) of every cubic
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
		if (bezierParams == null) {
			calculateDVals();
			calculatePsiVals();
			calculateThetaVals();
			calculatePhiVals();
			calculateCtrlPts();
		}
		return bezierParams;
	}

	/**
	 * Calculates the pairwise distances between the knots.
	 */
	private void calculateDVals() {
		// Skip last point if path is open
		int end = closed ? numPoints : numPoints - 1;
		for (int i = 0; i < end; i++) {
			Knot a = knots.get(i);
			Knot b = knots.get((i + 1) % numPoints);
			a.distance = a.c.dist(b.c);
			a.deltaX = b.c.x - a.c.x;
			a.deltaY = b.c.y - a.c.y;
		}
	}

	/**
	 * Calculates the turning angles ("psi") of the polyline which connect knots.
	 */
	private void calculatePsiVals() {
		// Skip first and last point if path is
		int start = closed ? 0 : 1;
		int end = closed ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			Knot prevKnot = knots.get(Math.floorMod(i - 1, numPoints));
			Knot knot = knots.get(i);
			double sin = prevKnot.deltaY / prevKnot.distance;
			double cos = prevKnot.deltaX / prevKnot.distance;

			knot.psi = Math.atan2(knot.deltaY * cos - knot.deltaX * sin, knot.deltaX * cos + knot.deltaY * sin);
		}
	}

	/**
	 * Creates five vectors which are coefficients of a linear system. Solving the
	 * system finds the value for theta (departure angle) at each point.
	 */
	private void calculateThetaVals() {
		/*
		 * Arrays holding the subdiagonals of the linear system that has to be solved to
		 * find the angles of the control points.
		 */
		final RealVector A = new ArrayRealVector(numPoints);
		final RealVector B = new ArrayRealVector(numPoints);
		final RealVector C = new ArrayRealVector(numPoints);
		final RealVector D = new ArrayRealVector(numPoints);
		final RealVector R = new ArrayRealVector(numPoints);

		// Calculate the entries of the five vectors.
		// Skip first and last point if path is open (no connecting bezier).
		int start = closed ? 0 : 1;
		int end = closed ? numPoints : numPoints - 1;
		for (int i = start; i < end; i++) {
			Knot z_h = (i == 0 ? knots.get(numPoints - 1) : knots.get(i - 1)); // prev
			Knot z_i = knots.get(i); // current
			Knot z_j = knots.get((i + 1) % numPoints); // next

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

		/*
		 * Special formulas for first and last rows of M (beziers having open paths),
		 * which don't follow the general rule.
		 */
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

		/*
		 * Solves a linear system to find departure angles (theta) at each knot. Since
		 * we already have turning angles each point (psi), the arrival angles (phi) can
		 * be obtained, since theta + phi + psi = 0 at each knot.
		 */
		DecompositionSolver solver = new LUDecomposition(M).getSolver();
		RealVector thetas = solver.solve(R);
		// Assign theta values to each Knot
		for (int i = 0; i < numPoints; i++) {
			knots.get(i).theta = thetas.getEntry(i);
		}

	}

	/**
	 * Calculates the arrival angles (phi) values for each knot via the
	 * relationship: <code>theta + phi + psi = 0</code>.
	 */
	private void calculatePhiVals() {
		for (Knot point : knots) {
			point.phi = -(point.psi + point.theta);
		}
	}

	/** Calculates the Bezier control points from knot z_i to z_{i+1}. */
	private void calculateCtrlPts() {
		bezierParams = new double[(closed ? numPoints : numPoints - 1)][8]; // Skip last point if path is open
		for (int i = 0; i < bezierParams.length; i++) {
			final Knot a = knots.get(i);
			final Knot b = knots.get((i + 1) % numPoints);
			final double rho = a.alpha * velocity(a.theta, b.phi); // coefficient
			final double sigma = b.beta * velocity(b.phi, a.theta); // coefficient

			Coordinate cp1 = a.c.add(b.c.sub(a.c).mult(rho / 3).mult(new Coordinate(Math.cos(a.theta), Math.sin(a.theta))));
			Coordinate cp2 = b.c.sub(b.c.sub(a.c).mult(sigma / 3).mult(new Coordinate(Math.cos(-b.phi), Math.sin(-b.phi))));

			int j = 0;
			bezierParams[i][j++] = a.c.x;
			bezierParams[i][j++] = a.c.y;
			bezierParams[i][j++] = cp1.x;
			bezierParams[i][j++] = cp1.y;
			bezierParams[i][j++] = cp2.x;
			bezierParams[i][j++] = cp2.y;
			bezierParams[i][j++] = b.c.x;
			bezierParams[i][j++] = b.c.y;
		}
	}

	/** Metafont's velocity function. */
	private static double velocity(final double theta, final double phi) {
		final double sinTheta = Math.sin(theta);
		final double sinPhi = Math.sin(phi);
		final double cosTheta = Math.cos(theta);
		final double cosPhi = Math.cos(phi);
		double numerator = 2 + Math.sqrt(2) * (sinTheta - sinPhi / 16) * (sinPhi - sinTheta / 16) * (cosTheta - cosPhi);
		double denominator = (1 + (Math.sqrt(5) - 1) / 2 * cosTheta + (3 - Math.sqrt(5)) / 2 * cosPhi);
		return numerator / denominator;
	}

}
