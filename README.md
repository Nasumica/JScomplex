# JScomplex
JavaScript complex arithmetic class with complex canvas exstension..

`var ef = new complex(Math.PI).muli.exp.add(1); //` Euler formula:  <img src="https://latex.codecogs.com/svg.latex?e^{i\pi}+1=0">

* Full set of complex arithmetic, elementary, hyperbolic and trigonometric functions and inverse.

* Complex polynomial equation solver. To solve<br>
z⁴ - (7 + 6 i) z³ - (1 - 30 i) z² + (67 - 4 i) z - (50 + 65 i) = 10 + 15 i<br>
enter code<br>
`
        var roots = []; new complex(10, 15).polysolve(roots, 1, [-7, -6], [-1, 30], [67, -4], [-50, -65]);
`<br>
and give roots: 4, 2 + i, i - 2, 3 + 4 i.

* Point manipulation operators: translation, rotation, projection, reflection;
isotomic, isogonal and harmonic conjugates;
directed distance;
line-line, line-circle and circle-circle intersection; circle tangent; bisection.

* Random complex numbers: bi-variate deviates for uniform, rectangular, polar, exponential, normal, Poisson distribution.

* Extension of CanvasRenderingContext2D prototype to using complex numbers as parameters.
Linkage call as: 
`ctx.begin.zmoveTo(z1).zlineTo(z2).zlineTo(z3).close.stroke();`

* Triangle centers solver using cartesian, trilinear, barycentric and tripolar coordinates.<br>
`var t = new complex(500, 300).triside(250, 350, 400, 15/180 * Math.PI);`<br>
creates triangle inclined by 15° with sides a = 250, b = 350, c = 400, and circumcenter at (500, 300).<br>
`ctx.begin.triangle(t).stroke();`<br>
draws him. To draw incircle enter code:<br>
`ctx.begin.ccircle(t.I).stroke();`<br>
or circumcircle<br>
`ctx.begin.ccircle(t.O).stroke();`<br>
or circumellipse<br>
`ctx.begin.eellipse(t.Oe).stroke();`<br>
Mark centroid G, incenter I, orthocenter H, circumcenter O and circumellipse foci F<sub>1</sub> and F<sub>2</sub><br>
`[t.G, t.I, t.H, t.O, t.Oe.F1, t.Oe.F2].map(function(z){ctx.begin.zcircle(z, 2).fill();});`

* Spiral of Theodorus, arithmetic and logarithmic spiral, Nautilus spiral, sunflower seeds spiral.
Implementation of superellipse. To draw Piet Hein supercircle (shape = 5/2) with center at (300, 200) and radius 150 enter<br>
`ctx.begin.supercircle(5/2, 300, 200, 150).stroke();`<br>
To draw 1000 sunflower seeds with center at (300, 200) enter something like this<br>
`for (var n = 1; n <= 1000; n++) ctx.begin.zcircle(new complex(300, 200).seed(n, 4), 2).fill();`

