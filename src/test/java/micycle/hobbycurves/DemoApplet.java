package micycle.hobbycurves;

import java.util.ArrayList;
import java.util.List;
import processing.core.PApplet;

public class DemoApplet extends PApplet {

	public static void main(String[] args) {
		PApplet.main(DemoApplet.class);
	}

	List<double[]> points;

	@Override
	public void settings() {
		size(800, 800);
	}

	@Override
	public void setup() {
		noFill();

		points = new ArrayList<>();
		for (int i = 0; i < 8; i++) {
			points.add(new double[] { random(100, width - 100), random(100, height - 100) });
		}
	}

	@Override
	public void draw() {
		background(255);

		double[] tensions = new double[points.size()];
		for (int i = 0; i < tensions.length; i++) {
			tensions[i] = random(1, 2);
		}

		IncrementalHobbyCurve ihc = new IncrementalHobbyCurve();
		points.forEach(p -> ihc.addKnot(p));

		if (points.size() >= 3) {
			boolean closed = false;
			float tension = map(mouseX, 0, width, 0.8f, 3);
			HobbyCurve hobbyCurve = new HobbyCurve(points.toArray(new double[0][0]), tension, closed, 0.25, 0.25);

			stroke(color(0, 165, 255, 100)); // blue
			strokeWeight(8);
			for (double[] b : ihc.getBeziers()) {
				int i = 0;
				bezier(b[i++], b[i++], b[i++], b[i++], b[i++], b[i++], b[i++], b[i++]);
			}

			stroke(color(255, 165, 0, 100)); // orange

			for (double[] b : hobbyCurve.getBeziers()) {
				int i = 0;
				bezier(b[i++], b[i++], b[i++], b[i++], b[i++], b[i++], b[i++], b[i++]);
			}
		}

		// processing curve
//		stroke(color(0, 50, 80, 200));
//		noFill();
//		beginShape();
//		curveVertex((float) points.get(0)[0], (float) points.get(0)[1]);
//		points.forEach(p -> curveVertex((float) p[0], (float) p[1]));
//		curveVertex((float) points.get(points.size() - 1)[0], (float) points.get(points.size() - 1)[1]);
//		endShape();

		stroke(color(255, 0, 0));
		strokeWeight(12);
		points.forEach(p -> point(p));

	}

	@Override
	public void mouseClicked() {
		points.add(new double[] { mouseX, mouseY });
	}

	@Override
	public void keyPressed() {
		points.clear();
		for (int i = 0; i < (int) random(5, 50); i++) {
			points.add(new double[] { random(100, width - 100), random(100, height - 100) });
		}
	}

	private void bezier(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
		super.bezier((float) x1, (float) y1, (float) x2, (float) y2, (float) x3, (float) y3, (float) x4, (float) y4);
	}

	private void point(double[] p) {
		super.point((float) p[0], (float) p[1]);
	}

}
