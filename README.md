# JScomplex
JavaScript complex arithmetic class with complex canvas exstension.

`var ef = new complex(Math.PI).muli.exp.add(1); //` Euler formula: e ^ (π * i) + 1

* Full set of complex arithmetic, elementary, hyperbolic and trigonometric functions.

* Basic linear algebra.

* Point manipulation operators; moving, rotation, projection, reflection;
isotonic, isogonal and harmonic conjugates;
directed distance;
line-line, line-circle and circle-circle intersection, circle tangent, bisection.

* Random complex numbers. Bi-variate deviates for uniform, polar, exponential, normal, Poisson distribution.

* Triangle centers solver using cartesian, trilinear, barycentric and tripolar coordinates.

* Arithmetic and logarithmic spiral, Nautilus spiral, sunflowers seed spiral.
Implementation of Johan Gielis superellipse (in 3 lines of code).

* Extension of CanvasRenderingContext2D prototype to using complex numbers as parameters.
Linkage call as: 

`ctx.startPath().zmoveTo(z1).zlineTo(z2).zlineTo(z3).endPath().stroke();`
