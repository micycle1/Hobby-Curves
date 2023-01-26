package micycle.hobbycurves;

import java.util.ArrayList;
import java.util.List;

/**
 * A variant of Hobby’s algorithm that builds up the curve one bezier at a time,
 * allowing for adding knots incrementally.
 * <p>
 * The incremental Hobby Curve is less “ideal” than the original (global)
 * variant, but has the property that adding new points will not affect earlier
 * segments.
 * 
 * @author Michael Carleton
 *
 */
public class IncrementalHobbyCurve implements IHobbyCurve {

	// implements https://github.com/loopspace/jsHobby/blob/master/hobby.js

	private static final double ha = Math.sqrt(2);
	private static final double hb = 1. / 16;
	private static final double hc = (3 - Math.sqrt(5)) / 2;
	private static final double hd = 1 - hc;

	private List<double[]> points;
	private List<double[]> beziers;

	private double thLast = Double.NaN;
	double[] lastBezier;

	/**
	 * Constructs a Incremental Hobby Curve having no initial knots.
	 */
	public IncrementalHobbyCurve() {
		points = new ArrayList<>();
		beziers = new ArrayList<>();
	}

	/**
	 * Constructs a Incremental Hobby Curve having initial knots.
	 */
	public IncrementalHobbyCurve(double[][] initialKnots) {
		this();
		for (double[] point : initialKnots) {
			addKnot(point);
		}
	}

	/**
	 * Adds a knot (k) this incremental Hobby Curve and returns the Bezier curve
	 * made by the previous two knots: k-2 and k-1.
	 * <p>
	 * Note: If the curve has less than 3 total knots, this method adds a knot but
	 * returns null, since 3 knots are required to create the first bezier curve.
	 * Furthermore, the bezier curve between k and k-1 is only accessible
	 * 
	 * @param point The point to be added and used to create the Bezier curve.
	 * 
	 * @return The Bezier curve created using the point, or null if there are less
	 *         than 3 points in the list.
	 */
	public double[] addKnot(double[] point) { // return bezier made by this point
		points.add(point);
		if (points.size() < 3) {
			return null;
		}
		HobbyResult h;
		if (points.size() == 3) {
			double[] pA = points.get(0);
			double[] pB = points.get(1);
			h = quickHobby(pA, pB, point, thLast);
		} else {
			double[] pA = points.get(points.size() - 3);
			double[] pB = points.get(points.size() - 2);
			h = quickHobby(pA, pB, point, thLast);
		}

		double[] bezier = flatten(h.a); // keep first curve only
		beziers.add(bezier);
		thLast = h.th;
		lastBezier = flatten(h.b);
		return bezier;
	}

	public double[][] getBeziers() {
		if (points.size() < 3) {
			return new double[0][0];
		}
		final double[][] bezierParams = new double[beziers.size() + 1][8];
		beziers.toArray(bezierParams);
		bezierParams[bezierParams.length - 1] = lastBezier; // append last bezier
		return bezierParams;
	}

	/**
	 * Flattens 2D [[ap1], [cp1]...] into 1D [ap1.x, ap1.y, cp1.x...]
	 * 
	 * @param b 2D array of bezier parameters
	 * @return flattened view of b
	 */
	private static double[] flatten(double[][] b) {
		return new double[] { b[0][0], b[0][1], b[1][0], b[1][1], b[2][0], b[2][1], b[3][0], b[3][1] };
	}

	/**
	 * Employ Hobby’s algorithm on 3 knots: a, b, c, to provide two cubic Bezier
	 * curves a->b and b->c. Of these, we later keep a->b and use that for the path
	 * between the kth and k+1st points.
	 * <p>
	 * We also remember the outgoing angle of the first segment and use that as the
	 * incoming angle on the next computation. The very first bezier of the hobby
	 * curve has no incoming angle (marked by being NaN).
	 */
	private static HobbyResult quickHobby(double[] a, double[] b, double[] c, double tha) {
		double da = dist(b, a);
		double db = dist(c, b);
		double wa = vecAng(b, a);
		double wb = vecAng(c, b);
		double psi = wb - wa;
		double thb, phb, phc;
		if (!Double.isNaN(tha)) {
			thb = -(2 * psi + tha) * db / (2 * db + da);
			phb = -psi - thb;
			phc = thb;
		} else {
			thb = -psi * db / (da + db);
			tha = -psi - thb;
			phb = tha;
			phc = thb;
		}
		return new HobbyResult(hobbySegment(a, tha, phb, b), hobbySegment(b, thb, phc, c), thb);
	}

	/**
	 * Generates a bezier curve connecting a and b.
	 * 
	 * @param a   knot a
	 * @param tha
	 * @param phb
	 * @param b   knot b
	 * @return
	 */
	private static double[][] hobbySegment(double[] a, double tha, double phb, double[] b) {
		double[] c = new double[] { b[0] - a[0], b[1] - a[1] };
		double sth = Math.sin(tha);
		double cth = Math.cos(tha);
		double sph = Math.sin(phb);
		double cph = Math.cos(phb);
		double alpha = ha * (sth - hb * sph) * (sph - hb * sth) * (cth - cph);
		double rho = (2 + alpha) / (1 + hd * cth + hc * cph);
		double sigma = (2 - alpha) / (1 + hd * cph + hc * cth);
		double[] ca = new double[] { a[0] + rho * (cth * c[0] - sth * c[1]) / 3, a[1] + rho * (sth * c[0] + cth * c[1]) / 3 };
		double[] cb = new double[] { b[0] - sigma * (cph * c[0] + sph * c[1]) / 3, b[1] - sigma * (-sph * c[0] + cph * c[1]) / 3 };
		return new double[][] { a, ca, cb, b };
	}

	private static double dist(double[] a, double[] b) {
		double dx = b[0] - a[0];
		double dy = b[1] - a[1];
		return Math.sqrt(dx * dx + dy * dy);
	}

	private static double vecAng(double[] a, double[] b) {
		double x = b[0] - a[0];
		double y = b[1] - a[1];
		return Math.atan2(y, x);
	}

	static class HobbyResult {

		double[][] a;
		double[][] b;
		double th;

		private HobbyResult(double[][] a, double[][] b, double th) {
			this.a = a;
			this.b = b;
			this.th = th;
		}

	}

}
