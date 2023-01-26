package micycle.hobbycurves;

/**
 * Models a point on the curve, having auxiliary parameters about the curve
 * passing through it.
 */
class Knot {

	/** Knot coordinates. */
	final Coordinate c;
	/** Tension of curve at this knot. */
	final double alpha;
	final double beta;
	/** Distance between this and next knots. */
	double distance = 0.0;
	double deltaX = 0.0;
	double deltaY = 0.0;
	/** Curve <b>departure</b> angle (radians) at this knot. */
	double theta = 0.0;
	/** Curve <b>arrival</b> angle (radians) at this knot. */
	double phi = 0.0;
	/** Curve <b>turning</b> angle (radians) at this knot. */
	double psi = 0.0;

	Knot(Coordinate complex, double alpha, double beta) {
		this.c = complex;
		this.alpha = alpha;
		this.beta = beta;
	}

	static class Coordinate {

		final double x, y;

		public Coordinate(double x, double y) {
			super();
			this.x = x;
			this.y = y;
		}

		double dist(Coordinate o) {
			double dx = x - o.x;
			double dy = y - o.y;
			return Math.sqrt(dx * dx + dy * dy);
		}

		double angle(Coordinate o) {
			double x = o.x - this.x;
			double y = o.y - this.y;
			return Math.atan2(y, x);
		}

		Coordinate mult(double factor) {
			return new Coordinate(x * factor, y * factor);
		}

		Coordinate mult(Coordinate o) {
			return new Coordinate(x * o.x - y * o.y, x * o.y + y * o.x);
		}

		Coordinate add(Coordinate o) {
			return new Coordinate(x + o.x, y + o.y);
		}

		Coordinate sub(Coordinate o) {
			return new Coordinate(x - o.x, y - o.y);
		}

	}

}