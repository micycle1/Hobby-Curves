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
import org.apache.commons.math3.util.Pair;

public class HobbyCurve {

	private List<HobbyPoint> points;
	List<Pair<Double, Double>> ctrl_pts;
	private final boolean cyclic;
	private final double begin_curl;
	private final double end_curl;
	private final int num_points;

	public HobbyCurve(List<Pair<Double, Double>> input_points, double tension, boolean cyclic, double begin_curl, double end_curl) {
		if (input_points.size() < 2) {
			throw new IllegalArgumentException("Hobby Algorithm needs more than two points");
		}
		this.points = new ArrayList<>(input_points.size());
		for (Pair<Double, Double> point : input_points) {
			HobbyPoint new_point = new HobbyPoint(new Complex(point.getFirst(), point.getSecond()), 1.0 / tension, 1.0 / tension);
			points.add(new_point);
		}
		this.cyclic = cyclic;
		this.begin_curl = begin_curl;
		this.end_curl = end_curl;
		this.num_points = points.size();
	}

	List<Pair<Double, Double>> getCtrlPts() {
		calculateDVals();
		calculatePsiVals();
		calculateThetaVals();
		calculatePhiVals();
		calculateCtrlPts();
		return ctrl_pts;
	}

	/** Calculates the pairwise distances between the points. */
	void calculateDVals() {
		// Skip last point if path is non-cyclic
		int end = cyclic ? num_points : num_points - 1;
		for (int i = 0; i < end; i++) {
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % num_points);
			z_i.d_val = z_i.cmplx.subtract(z_j.cmplx).abs();
		}
	}

	/** Calculates the psi values by subtracting pairwise phases. */
	void calculatePsiVals() {
		// Skip first and last point if path is non-cyclic
		int start = cyclic ? 0 : 1;
		int end = cyclic ? num_points : num_points - 1;
		for (int i = start; i < end; i++) {
			HobbyPoint z_h = (i == 0 ? points.get(num_points - 1) : points.get(i - 1));
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % num_points);
			Complex polygonal_turn = z_j.cmplx.subtract(z_i.cmplx).divide(z_i.cmplx.subtract(z_h.cmplx));
			z_i.psi = Math.atan2(polygonal_turn.getImaginary(), polygonal_turn.getReal());
		}
	}

	/**
	 * Calculates the theta values by creating a linear system whose solutions are
	 * the values.
	 */
	void calculateThetaVals() {
		final RealVector A = new ArrayRealVector(num_points);
		final RealVector B = new ArrayRealVector(num_points);
		final RealVector C = new ArrayRealVector(num_points);
		final RealVector D = new ArrayRealVector(num_points);
		final RealVector R = new ArrayRealVector(num_points);

		// Calculate the entries of the five vectors.
		// Skip first and last point if path is non-cyclic.
		int start = cyclic ? 0 : 1;
		int end = cyclic ? num_points : num_points - 1;
		for (int i = start; i < end; i++) {
			HobbyPoint z_h = (i == 0 ? points.get(num_points - 1) : points.get(i - 1));
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % num_points);

			A.setEntry(i, z_h.alpha / (Math.pow(z_i.beta, 2) * z_h.d_val));
			B.setEntry(i, (3 - z_h.alpha) / (Math.pow(z_i.beta, 2) * z_h.d_val));
			C.setEntry(i, (3 - z_j.beta) / (Math.pow(z_i.beta, 2) * z_h.d_val));
			D.setEntry(i, z_j.beta / (Math.pow(z_i.alpha, 2) * z_i.d_val));
			R.setEntry(i, -B.getEntry(i) * z_i.psi - D.getEntry(i) * z_j.psi);
		}

		RealMatrix M = new Array2DRowRealMatrix(num_points, num_points);
		// Set up matrix M such that the soln. Mx = R are the theta values.
		for (int i = start; i < end; i++) {
			// Fill i-th row of M
			if (i == 0) {
				M.setEntry(i, num_points - 1, A.getEntry(i));
			} else {
				M.setEntry(i, i - 1, A.getEntry(i));
			}
			M.setEntry(i, i, B.getEntry(i) + C.getEntry(i));
			M.setEntry(i, (i + 1) % num_points, D.getEntry(i));
		}

		// Special formulas for first and last rows of M with non-cyclic paths.
		if (!cyclic) {
			// First row of M
			double alpha_0 = points.get(0).alpha;
			double beta_1 = points.get(1).beta;
			double xi_0 = (Math.pow(alpha_0, 2) * begin_curl) / Math.pow(beta_1, 2);
			M.setEntry(0, 0, alpha_0 * xi_0 + 3 - beta_1);
			M.setEntry(0, 1, (3 - alpha_0) * xi_0 + beta_1);
			R.setEntry(0, -((3 - alpha_0) * xi_0 + beta_1) * points.get(1).psi);
			// Last row of M
			double alpha_n_1 = points.get(num_points - 2).alpha;
			double beta_n = points.get(num_points - 1).beta;
			double xi_n = (Math.pow(beta_n, 2) * end_curl) / Math.pow(alpha_n_1, 2);
			M.setEntry(num_points - 1, num_points - 2, (3 - beta_n) * xi_n + alpha_n_1);
			M.setEntry(num_points - 1, num_points - 1, (beta_n * xi_n + 3 - alpha_n_1));
			R.setEntry(num_points - 1, 0);
		}

		// Solve for theta values
		DecompositionSolver solver = new LUDecomposition(M).getSolver();
		RealVector thetas = solver.solve(R);
		// Assign theta values to each HobbyPoint
		for (int i = 0; i < num_points; i++) {
			points.get(i).theta = thetas.getEntry(i);
		}

	}

	/**
	 * Calculates the phi_k values via the relationship theta_k + phi_k + psi_k = 0.
	 */
	void calculatePhiVals() {
		for (HobbyPoint point : points) {
			point.phi = -(point.psi + point.theta);
		}
	}

	/** Calculates the Bezier control points from z_i to z_{i+1}. */
	void calculateCtrlPts() {
		ctrl_pts = new ArrayList<>(num_points);
		// Skip last point if path is non-cyclic
		int end = cyclic ? num_points : num_points - 1;
		for (int i = 0; i < end; i++) {
			HobbyPoint z_i = points.get(i);
			HobbyPoint z_j = points.get((i + 1) % num_points);
			double rho_coefficient = z_i.alpha * velocity(z_i.theta, z_j.phi);
			double sigma_coefficient = z_j.beta * velocity(z_j.phi, z_i.theta);

			Complex ctrl_point_a = z_i.cmplx.add(z_j.cmplx.subtract(z_i.cmplx).multiply(rho_coefficient / 3)
					.multiply(new Complex(Math.cos(z_i.theta), Math.sin(z_i.theta))));
			Complex ctrl_point_b = z_j.cmplx.subtract(z_j.cmplx.subtract(z_i.cmplx).multiply(sigma_coefficient / 3)
					.multiply(new Complex(Math.cos(-z_j.phi), Math.sin(-z_j.phi))));
			Pair<Double, Double> p = Pair.create(ctrl_point_a.getReal(), ctrl_point_a.getImaginary());
			ctrl_pts.add(p);
			ctrl_pts.add(Pair.create(ctrl_point_b.getReal(), ctrl_point_b.getImaginary()));
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
