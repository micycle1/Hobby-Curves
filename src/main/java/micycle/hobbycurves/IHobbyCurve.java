package micycle.hobbycurves;

public interface IHobbyCurve {

	/**
	 * Returns a list of parameters for the cubic bezier curves comprising this
	 * Hobby Curve.
	 * <p>
	 * Each curve has 4 parameters (each an (x,y) coordinate pair), so each row of
	 * the output has the form...
	 * <p>
	 * <code>[ap1.x, ap1.y, cp1.x ,cp1.y, cp2.x, cp2.y, ap2.x, ap2.y]</code>
	 * <p>
	 * ...where <code>ap</code> denotes "anchor point" and <code>cp</code> denotes
	 * "control point".
	 * 
	 * @return a list of bezier curve parameters
	 */
	double[][] getBeziers();

}
