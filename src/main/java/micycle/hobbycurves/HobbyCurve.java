package micycle.hobbycurves;

import java.util.ArrayList;
import java.util.List;

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
 *
 */
public class HobbyCurve implements IHobbyCurve {

	// based on https://github.com/ltrujello/Hobby_Curve_Algorithm/

	private List<Knot> knots;
	private double[][] bezierParams;
	private final boolean closed;
	private final double beginCurl;
	private final double endCurl;
	private final int numPoints;

	/**
	 * Constructor for creating a new "Hobby Curve".
	 *
	 * @param inputKnots a 2D array of [x, y] coordinates of the data knots that the
	 *                   curve should pass through.
	 * @param tension    a parameter that controls the tension of the curve's
	 *                   "knots". A value of 1 is a good starting point.
	 * @param closed     a boolean value indicating whether the curve should form a
	 *                   closed loop, joining the first and last knots
	 * @param beginCurl  a value that controls the amount of curl at the start of
	 *                   the curve (the angle of incoming bezier originating at the
	 *                   starting point). Accepts any value, but should generally be
	 *                   positive; a value of 1 is a good starting point. This value
	 *                   is only effective if the curve is unclosed.
	 * @param endCurl    a value that controls the amount of curl at the end of the
	 *                   curve (the angle out outgoing bezier terminating at the
	 *                   final point). Accepts any value, but should generally be
	 *                   positive; a value of 1 is a good starting point. This value
	 *                   is only effective if the curve is unclosed.
	 */
	public HobbyCurve(double[][] inputKnots, double tension, boolean closed, double beginCurl, double endCurl) {
		if (inputKnots.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two knots");
		}
		this.knots = new ArrayList<>(inputKnots.length);
		for (double[] point : inputKnots) {
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
	 * @param inputKnots a 2D array of [x, y] coordinates of the data knots that the
	 *                   curve should pass through.
	 * @param tensions   a list of tension parameters that control the tension of
	 *                   the curve at each corresponding knots. A value of 1 is a
	 *                   good starting point.
	 * @param closed     a boolean value indicating whether the curve should form a
	 *                   closed loop, joining the first and last knots
	 * @param beginCurl  a value that controls the amount of curl at the start of
	 *                   the curve (the angle of incoming bezier originating at the
	 *                   starting point). Accepts any value, but should generally be
	 *                   positive; a value of 1 is a good starting point. This value
	 *                   is only effective if the curve is unclosed.
	 * @param endCurl    a value that controls the amount of curl at the end of the
	 *                   curve (the angle out outgoing bezier terminating at the
	 *                   final point). Accepts any value, but should generally be
	 *                   positive; a value of 1 is a good starting point. This value
	 *                   is only effective if the curve is unclosed.
	 */
	public HobbyCurve(double[][] inputKnots, double[] tensions, boolean closed, double beginCurl, double endCurl) {
		if (inputKnots.length < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two knots");
		}
		if (inputKnots.length != tensions.length) {
			throw new IllegalArgumentException("Points and tensions arrays are different sizes. They should be equal.");
		}
		this.knots = new ArrayList<>(inputKnots.length);
		for (int i = 0; i < inputKnots.length; i++) {
			double[] point = inputKnots[i];
			double tension = tensions[i];
			Knot knot = new Knot(new Coordinate(point[0], point[1]), 1.0 / tension, 1.0 / tension);
			knots.add(knot);
		}
		this.closed = closed;
		this.beginCurl = beginCurl;
		this.endCurl = endCurl;
		this.numPoints = knots.size();
	}

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

	/** Original Hobby algorithm, but with O(n) linear solver. */
	private void calculateThetaVals() {

		final int n = numPoints;
		final boolean cyc = closed;

		/* --- build the three diagonals and the right–hand side ------------- */

		double[] low = new double[n]; // a_{i, i-1}
		double[] diag = new double[n]; // a_{i, i}
		double[] up = new double[n]; // a_{i, i+1}
		double[] rhs = new double[n];

		int start = cyc ? 0 : 1;
		int end = cyc ? n : n - 1;

		for (int i = start; i < end; i++) {
			Knot zh = knots.get(Math.floorMod(i - 1, n));
			Knot zi = knots.get(i);
			Knot zj = knots.get((i + 1) % n);

			double A = zh.alpha / (zi.beta * zi.beta * zh.distance);
			double B = (3 - zh.alpha) / (zi.beta * zi.beta * zh.distance);
			double C = (3 - zj.beta) / (zi.alpha * zi.alpha * zi.distance);
			double D = zj.beta / (zi.alpha * zi.alpha * zi.distance);

			low[i] = A;
			diag[i] = B + C;
			up[i] = D;
			rhs[i] = -B * zi.psi - D * zj.psi;
		}

		/* --- special first / last rows for the open case -------------------- */
		if (!cyc) {
			Knot z0 = knots.get(0), z1 = knots.get(1);
			double xi0 = (z0.alpha * z0.alpha * beginCurl) / (z1.beta * z1.beta);

			diag[0] = z0.alpha * xi0 + 3 - z1.beta;
			up[0] = (3 - z0.alpha) * xi0 + z1.beta;
			rhs[0] = -up[0] * z1.psi;

			Knot zn_1 = knots.get(n - 2), zn = knots.get(n - 1);
			double xin = (zn.beta * zn.beta * endCurl) / (zn_1.alpha * zn_1.alpha);

			low[n - 1] = (3 - zn.beta) * xin + zn_1.alpha;
			diag[n - 1] = zn.beta * xin + 3 - zn_1.alpha;
			// rhs[n-1] is already 0
		}

		// solve
		double[] theta = cyc ? solveCyclic(low, diag, up, rhs) : solveTri(low, diag, up, rhs);

		for (int i = 0; i < n; i++) {
			knots.get(i).theta = theta[i];
		}
	}

	/**
	 * /** Solve a×x = d where 'a' is tridiagonal (open curve).
	 */
	private static double[] solveTri(double[] low, double[] diag, double[] up, double[] rhs) {

		int n = diag.length;

		// Forward sweep -------------------------------------------------
		for (int i = 1; i < n; i++) {
			double m = low[i] / diag[i - 1];
			diag[i] -= m * up[i - 1];
			rhs[i] -= m * rhs[i - 1];
		}

		// Back substitution --------------------------------------------
		double[] x = new double[n];
		x[n - 1] = rhs[n - 1] / diag[n - 1];

		for (int i = n - 2; i >= 0; i--) {
			x[i] = (rhs[i] - up[i] * x[i + 1]) / diag[i];
		}
		return x;
	}

	/** Solve cyclic tridiagonal system (closed curve). */
	private static double[] solveCyclic(double[] low, double[] diag, double[] up, double[] rhs) {

		int n = diag.length;
		double[] u = new double[n];
		double[] v = new double[n];

		// u = [gamma, 0, 0, …, beta]
		double gamma = -diag[0]; // any non–zero number is fine
		u[0] = gamma;
		u[n - 1] = low[0];

		v[0] = 1.0;
		v[n - 1] = up[n - 1] / gamma;

		// Modify first and last diagonal entries
		double[] diag2 = diag.clone();
		diag2[0] -= gamma;
		diag2[n - 1] -= up[n - 1] * low[0] / gamma;

		// Solve two ordinary tridiagonal systems
		double[] y = solveTri(low, diag2, up, rhs.clone());
		double[] z = solveTri(low, diag2, up, u);

		// Combine with Sherman–Morrison
		double fact = (y[0] + up[n - 1] * y[n - 1] / gamma) / (1.0 + v[0] + up[n - 1] * z[n - 1] / gamma);

		double[] x = new double[n];
		for (int i = 0; i < n; i++) {
			x[i] = y[i] - fact * z[i];
		}
		return x;
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
