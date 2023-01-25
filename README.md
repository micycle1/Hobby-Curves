[![](https://jitpack.io/v/micycle1/Hobby-Curves.svg)](https://jitpack.io/#micycle1/Hobby-Curves)

# Hobby-Curves

_Hobby Curves_, in Java.

John Hobby’s algorithm [[1]] produces a smooth curve through a given set of points. The curve comprises a chain of cubic Bézier curves whose endpoints pass through the points. The parameters of the Bézier curves are chosen such that they join smoothly, forming one long curve through the point set.

Hobby Curves are more visually pleasing than curves produced via most other techniques.

This Java implementation is based on _Luke Trujillo's_ C++ [implementation](https://github.com/ltrujello/Hobby_Curve_Algorithm]).

## Example

```
double[][] points = new double[4][2];
points[0] = new double[] { 0.0, 0.0 };
points[1] = new double[] { 10.0, 15.0 };
points[2] = new double[] { 20.0, 0.0 };
points[3] = new double[] { 10.0, -10.0 };

double tension = 1;
double closed = true;
double beginCurl = 1;
double endCurl = 1;

HobbyCurve curve = new HobbyCurve(points, tension, closed, beginCurl, endCurl);

for (double[] bezier : curve.getBeziers()) {
	System.out.println(Arrays.toString(bezier));
}
```

https://mirror.apps.cam.ac.uk/pub/tex-archive/graphics/pgf/contrib/hobby/hobby_code.pdf
http://ctan.math.utah.edu/ctan/tex-archive/graphics/pgf/contrib/hobby/hobby.pdf

[1]: https://www.researchgate.net/publication/226514776_Smooth_easy_to_compute_interpolating_splines
