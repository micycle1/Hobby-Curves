# Hobby-Curves
John Hobby's curve drawing algorithm.

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
