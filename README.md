# JScomplex
JavaScript complex arithmetic class with complex canvas exstension..

`var ef = new complex(Math.PI).muli.exp.add(1); //` Euler formula:  <img src="https://latex.codecogs.com/svg.latex?e^{i\pi}+1=0">

* Full set of complex arithmetic, elementary, hyperbolic and trigonometric functions and inverse.

* Complex polynomial equation solver. To solve<br>
z⁴ - (7 + 6 i) z³ - (1 - 30 i) z² + (67 - 4 i) z - (50 + 65 i) = 10 + 15 i<br>
enter code<br>
`
        var z = new complex(10, 15).polysolve(1, [-7, -6], [-1, 30], [67, -4], [-50, -65]);
`<br>
and give roots: 4, 2 + i, i - 2, 3 + 4 i.

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
