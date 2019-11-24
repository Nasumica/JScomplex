# JScomplex
JavaScript complex arithmetic class with complex canvas exstension..

`var ef = new complex(Math.PI).muli.exp.add(1); //` Euler formula:  <img src="https://latex.codecogs.com/svg.latex?e^{i\pi}+1">

* Full set of complex arithmetic, elementary, hyperbolic and trigonometric functions and inverse.

* Point manipulation operators: translation, rotation, projection, reflection;
isotonic, isogonal and harmonic conjugates;
directed distance;
line-line, line-circle and circle-circle intersection; circle tangent; bisection.

* Random complex numbers: bi-variate deviates for uniform, rectangular, polar, exponential, normal, Poisson distribution.

* Triangle centers solver using cartesian, trilinear, barycentric and tripolar coordinates.

* Spiral of Theodorus, arithmetic and logarithmic spiral, Nautilus spiral, sunflower seeds spiral.
Implementation of superellipse.

* Extension of CanvasRenderingContext2D prototype to using complex numbers as parameters.
Linkage call as: 
`ctx.begin.zmoveTo(z1).zlineTo(z2).zlineTo(z3).close.stroke();`
