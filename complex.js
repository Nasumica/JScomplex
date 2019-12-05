class complex {
/*
	Complex arithmetic with canvas extension
	
	JavaScript implementation by Србислав Д. Нешић, November 2019
	srbislav.nesic@gmail.com
	Nasumica Agencija
	
	https://github.com/Nasumica/JScomplex

	Euler identity e ^ (i * π) + 1 = 0: new complex(Math.PI).muli.exp.add(1);
*/
	constructor(u, v){
	/*
		var z = new complex();                   // z = 0
		var z = new complex(a);                  // z = a
		var z = new complex(a, b);               // z = a + i b
		var z = new complex({});                 // z = 0
		var z = new complex({x: 3});             // z = 3
		var z = new complex({y: 4});             // z = 4 i
		var z = new complex({x: 3, y: 4});       // z = 3 + 4 i
		var z = new complex({r: 2, y: 4, x: 3}); // z = 3 + 4 i (r is ignored)
		var z = new complex([]);                 // z = 0
		var z = new complex([5]);                // z = 5
		var z = new complex([6, 7]);             // z = 6 + 7 i
		var z = new complex([{x: 1, y: 2}]);     // z = NaN
	*/
		this.x = 0; this.y = 0;
		if (arguments.length > 1) this.xiy(u, v); else
		if (arguments.length > 0) 
		if (typeof u === 'object') {
			if (Array.isArray(u)){
				if (u.length > 1) this.xiy(u[0], u[1]); else
				if (u.length > 0) this.xiy(u[0]);
			} else this.asg(u); 
		} else this.xiy(u);
	}
	nop(){// no operation
		return this;
	}
	get z(){// as object
		return {x: this.x, y: this.y};
	}
	xiy(x = 0, y = 0){// x + i y
		this.x = Number(x); this.y = Number(y); return this;
	}
	asg(z, condition = true){// copy z to this
		if (condition){
			this.zero;
			if (typeof z.x !== 'undefined') this.x = Number(z.x);
			if (typeof z.y !== 'undefined') this.y = Number(z.y);
		}
		return this;
	}
	cis(rho, theta){// ρ cis θ = ρ (cos θ + i sin θ) = ρ e^(iθ)
		if (arguments.length < 2) {theta = rho; rho = 1;} // if ρ is ommited then ρ = 1
		return this.asg(cis(theta, rho));
	}
	iff(yes, no, condition = true){
		return this.asg(condition ? yes : no);
	}
	obj(z){// copy this to z
		z.x = this.x; z.y = this.y; return this;
	}
	swap(z){// swap values with z
		var t = {x: z.x, y: z.y}; return this.obj(z).asg(t);
	}
	get sqrabs(){// |z|² = z * z' = ρ²
		return this.x * this.x + this.y * this.y;
	}
	get abs(){// magnitude ρ
		return this.x == 0 || this.y == 0 ? abs(this.x + this.y) : sqrt(this.sqrabs);
	}
	get arg(){// angle θ
		return atan(this.y, this.x);
	}
	isEq(x, y = 0){
		return this.x == x && this.y == y;
	}
	is(z){
		return this.isEq(z.x, z.y);
	}
	get is0(){
		return this.isEq(0);
	}
	get is1(){
		return this.isEq(1);
	}
	get isI(){
		return this.isEq(0, 1);
	}
	get isNum(){
		return Number.isFinite(this.x) && Number.isFinite(this.y);
	}
	get isInf(){
		return abs(this.x) === Infinity || abs(this.y) === Infinity;
	}
	get isNaN(){
		return Number.isNaN(this.x) || Number.isNaN(this.y);
	}
	isEps(z = {x: 0, y: 0}, eps = 1e-15){// with respect to epsilon
		return abs(this.x - z.x) + abs(this.y - z.y) < eps; // must improve
	}
	get zero(){
		return this.xiy(0);
	}
	get one(){
		return this.xiy(1);
	}
	get i(){// sqrt(-1)
		return this.xiy(0, 1);
	}
	get inf(){ // complex infinity
		return this.xiy(1/0, 1/0);
	}
	get nan(){ // not a number
		return this.xiy(0/0, 0/0);
	}
	get real(){// projection to x-axis
		return this.xiy(this.x, 0);
	}
	get imag(){// projection to y-axis
		return this.xiy(0, this.y);
	}
	get conjg(){// reflection to x-axis
		return this.xiy(this.x, -this.y);
	}
	get pos(){// 1st quadrant
		return this.xiy(abs(this.x), abs(this.y));
	}
	get neg(){// -z = z rotate 180°
		return this.xiy(-this.x, -this.y);
	}
	get muli(){// z * i = z rotate 90°
		return this.xiy(-this.y, this.x);
	}
	get divi(){// z / i = z rotate -90° = - i * z
		return this.xiy(this.y, -this.x);
	}
	get inv(){// reflection to line x = y
		return this.xiy(this.y, this.x);
	}
	rect(z){// rectangle area
		return (this.x - z.x) * (this.y - z.y);
	}
	trap(z){// signed trapezoid area
		return (this.x - z.x) * (this.y + z.y) / 2;
	}
	adv(x, y = x){
		return this.xiy(this.x + x, this.y + y);
	}
	zadv(z){// advance
		return this.adv(z.x, z.y);
	}
	scl(x, y = x){
		return this.xiy(this.x * x, this.y * y);
	}
	zscl(z){// scale up
		return this.scl(z.x, z.y);
	}
	lcs(x, y = x){
		return this.xiy(this.x / x, this.y / y);
	}
	zlcs(z){// scale down
		return this.lsc(z.x, z.y);
	}
	add(x, y = 0){
		return this.xiy(this.x + x, this.y + y);
	}
	zadd(z){// this + z
		return this.add(z.x, z.y);
	}
	sub(x, y = 0){
		return this.xiy(this.x - x, this.y - y);
	}
	zsub(z){// this - z
		return this.sub(z.x, z.y);
	}
	bus(x, y = 0){
		return this.sub(x, y).neg;
	}
	zbus(z){// z - this
		return this.bus(z.x, z.y);
	}
	mul(x, y = 0){
		return this.xiy(this.x * x - this.y * y, this.x * y + this.y * x);
	}
	zmul(z){ // this * z
		return this.mul(z.x, z.y);
	}
	div(x, y = 0){
		if (y == 0) return this.lcs(x); else 
		if (x == 0) return this.lcs(y).divi; else
			return this.mul(x, -y).lcs(x * x + y * y);
	}
	zdiv(z){// this / z
		return this.div(z.x, z.y);
	}
	vid(x, y = 0){
		return this.asg(new complex(x, y).zdiv(this));
	}
	zvid(z){// z / this
		return this.vid(z.x, z.y);
	}
	get recip(){// 1 / this
		if (this.y == 0) return this.xiy(1 / this.x); else 
		if (this.x == 0) return this.xiy(1 / this.y).divi; else
			return this.conjg.lcs(this.sqrabs);
	}
	mod(x = 1, y = 0){
		var z = new complex(this);
		return this.div(x, y).floor().mul(x, y).zbus(z);
	}
	zmod(z){
		return this.mod(z.x, z.y);
	}
	get sqr(){// this²
		return this.zmul(this);
	}
	get cub(){// this³
		return this.zmul(new complex(this).sqr);
	}
	get sqrt(){
		if (this.is0 || this.is1) return this; else
			if (this.y == 0){
				if (this.x < 0)
					return this.xiy(sqrt(-this.x)).muli;
				else
					return this.xiy(sqrt(this.x));
			} else return this.cis(sqrt(this.abs), this.arg/2);
	}
	get cbrt(){
		if (this.is0 || this.is1) return this; else
			if (this.y == 0)
				return this.xiy(Math.sign(this.x) * pow(abs(this.x), 1/3));
			else
				return this.root(3);
	}
	get unit(){// unit vector
		var d = this.abs; return d == 0 ? this : this.div(d);
	}
	get theo(){// spiral of Theodorus
		// a = 2*sqrt(n) - 2.1577829966594462209291427868295777235
		// f(n + 1) = (1 + i / sqrt(n + 1)) * f(n)
		if (this.isI) this.xiy(0); else
		if (this.is0) this.xiy(1); else
			this.zadd(new complex(this).unit.muli);
		return this;
	}
	inc(v = 1){// ρ = |ρ + v|
		return this.is0 ? this.xiy(v) : this.zadd(new complex(this).unit.mul(v));
	}
	dec(v = 1){// ρ = |ρ - v|
		return this.inc(-v);
	}
	round(r = 1){
		return this.scl(r).xiy(Math.round(this.x), Math.round(this.y)).lcs(r); 
	}
	floor(r = 1){
		return this.scl(r).xiy(Math.floor(this.x), Math.floor(this.y)).lcs(r); 
	}
	ceil(r = 1){
		return this.scl(r).xiy(Math.ceil(this.x), Math.ceil(this.y)).lcs(r); 
	}
	get integer(){
		return this.round();
	}
	dcp(x, y = 0){// Re = dot product, Im = cross product (not commutative)
		return this.conjg.mul(x, y);
	}
	zdcp(z){
		return this.dcp(z.x, z.y);
	}
	get exp(){// e^(x + iy) = e^x * e^(iy) = e^x * (cos y + i sin y) = e^x cis y
		if (this.isEq(0, pi)) return this.xiy(-1); else // to Euler
		return this.cis(exp(this.x), this.y);
	}
	get log(){// ln(ρ cis θ) = ln(ρ * e^(iθ)) = ln ρ + ln(e^(iθ)) = ln ρ + iθ
		if (this.isEq(-1)) return this.xiy(0, pi); else // ditto
		return this.xiy(log(this.abs), this.arg);
	}
	pow(x, y = 0){
		if (this.is1) return this; else
		if (y == 0 && x ==  0) return this.one; else
		if (y == 0 && x ==  1) return this; else
		if (y == 0 && x == -1) return this.recip; else
		if (y == 0 && x > 0 && this.is0) return this; else
			return this.log.mul(x, y).exp;
	}
	zpow(z){
		return this.pow(z.x, z.y);
	}
	npow(n){// thisⁿ (integer n); De Moivre's formula: (cis θ)ⁿ = cis (θ·n)
		var z = new complex(this);
		var i = Math.floor(abs(n));
		this.one;
		while (i > 0){
			if (i % 2 == 1) this.zmul(z);
			z.sqr;  i >>>= 1;
		}
		if (n < 0) this.recip;
		return this;
	}
	root(x, y = 0){
		if (this.is1) return this; else
		if (y == 0 && x ==  1) return this; else
		if (y == 0 && x == -1) return this.recip; else
		if (y == 0 && x > 0 && this.is0) return this; else
			return this.log.div(x, y).exp;
	}
	zroot(z){
		return this.root(z.x, z.y);
	}
	get sinh(){// (e^z - e^-z)/2
		this.exp; return this.zsub(new complex(this).recip).div(2);
	}
	get cosh(){// (e^z + e^-z)/2
		this.exp; return this.zadd(new complex(this).recip).div(2);
	}
	get tanh(){// (e^2z - 1)/(e^2z + 1)
		this.mul(2).exp; var z = new complex(this).add(1);
		return this.sub(1).zdiv(z);
	}
	get sin(){// i * sinh(z) = sin(i * z) => sin(z) = sinh(z / i) * i
		return this.divi.sinh.muli;
	}
	get cos(){// cosh(z) = cos(i * z) => cos(z) = cosh(z / i)
		return this.divi.cosh;
	}
	get tan(){// i * tanh(z) = tan(i * z)
		return this.divi.tanh.muli;
	}
	get sinc(){// sin(z)/z
		var z = new complex(this);
		return this.is0 ? this.one : this.sin.zdiv(z);
	}
	get asinh(){// ln(z + sqrt(z² + 1))
		var z = new complex(this);
		return this.sqr.add(1).sqrt.zadd(z).log;
	}
	get acosh(){// ln(z + sqrt(z² - 1))
		var z = new complex(this);
		return this.sqr.sub(1).sqrt.zadd(z).log;
	}
	get atanh(){// ln((1 + z)/(1 - z))/2 = (ln(1 + z) - ln(1 - z))/2
		var z = new complex(1).zsub(this);
		return this.add(1).zdiv(z).log.div(2);
	}
	get asin(){
		return this.muli.asinh.divi;
	}
	get acos(){
		return this.acosh.divi;
	}
	get atan(){
		return this.muli.atanh.divi;
	}
	horner(){// Horner's scheme polynom evaluate
		var z = {}; this.obj(z).zero;
		for (var i = 0; i < arguments.length; i++)
			this.zmul(z).zadd(new complex(arguments[i]));
		return this;
	}
	spline(){// first arg for this = 0, last arg for this = 1, else somewhere between smoothly
		if (arguments.length > 0){
			var t = new complex(this), s = new complex(1).zsub(t); // t = this, s = 1 - t
			var n = arguments.hi, m = 0, b = 1;  this.zero;
			while (n >= 0){// $b = {arguments.hi \choose n}$
				this.zadd(new complex(arguments[n]).mul(b).zmul(new complex(s).npow(m)).zmul(new complex(t).npow(n)));
				m++; b *= n; b /= m; n--;
			}
		}
		return this;
	}
	get smoothstep(){// smoothstep real part
		var x = this.x;
		if (x >=  1) this.x = 1; else
		if (x <=  0) this.x = 0; else
		if (x > 0.5) this.bus(1).smoothstep.bus(1); else
		if (x < 0.5) this.x = x * x * (3 - 2 * x); // Hermite cubic spline
		return this;
	}
	rotate(angle){// rotate this about origin by given angle
		var h = pi/2; // few exact values
		if (angle == 0) return this; else 
		if (angle == h) return this.muli; else 
		if (angle == -h || angle == 3 * h) return this.divi; else 
		if (abs(angle) == pi) return this.neg; else
			return this.zmul(cis(angle));
	}
	zrotate(z){// rotate this about origin by angle of vector (0, 0)--z
		var h = z.x * z.x + z.y * z.y;
		return h == 0 ? this : this.zmul(z).div(sqrt(h));
	}
	vrotate(z1, z2){// rotate this about origin by angle of vector z1--z2
		return this.zrotate(new complex(z2).zsub(z1));
	}
	about(x, y, angle){// rotate this about point (x, y)
		return this.sub(x, y).rotate(angle).add(x, y);
	}
	zabout(z, angle){// rotate this about point z
		return this.about(z.x, z.y, angle);
	}
	cyclic(m, n = 1){// rotate by m/n part of full circle
		return this.rotate(2 * pi * m/n);
	}
	times(z1, z2 = {x: 0, y: 0}){// this = z1 + times * (z2 - z1)
		return this.zsub(z1).zdiv(new complex(z2).zsub(z1)); 
	}
	inter(z1, z2 = {x: 0, y: 0}){// inversion of times; t = this, s = 1-t, result = z1*s + z2*t
		//return this.spline(z1, z2);
		return this.zmul(new complex(z2).zsub(z1)).zadd(z1);
	}
	quadraticeq(A, B, C){// quadratic equation solver A z² + B z + C = this
		this.z1 = {}; this.z2 = {}; // result
		var a = new complex(A); // a = A
		var b = new complex(B); // b = B
		var c = new complex(C).zsub(this); // c = C - this
		if (a.is0)// B z + C = 0
			c.zdiv(b).neg.obj(this.z1).obj(this.z2); // z1 = z2 = -C/B
		else if (b.is0)// A z² + C = 0
			c.zdiv(a).neg.sqrt.obj(this.z1).neg.obj(this.z2); // -z2 = z1 = sqrt(-C/A)
		else if (c.is0)// (A z + B) z = 0
			b.zdiv(a).neg.obj(this.z1).zero.obj(this.z2); // z1 = -B/A, z2 = 0
		else {// A z² + B z + C = 0
			var d = new complex(b.neg).sqr.zsub(c.mul(2).zmul(a.mul(2))).sqrt;
			// b = -B; a = 2 A; c = 4 A C; d = sqrt(B² - 4 A C); 
			c.asg(b).zadd(d).zdiv(a).obj(this.z1);
			c.asg(b).zsub(d).zdiv(a).obj(this.z2);
		}
		return this;
	}
	quadtimes(z1, z2, z3){// similar to linear times but quadratic
		// solve (z1 + z3 - 2 z2) t² + 2 (z2 - z1) t + z1 = this
		this.quadraticeq(
			new complex(z2).mul(-2).zadd(z1).zadd(z3), // A = z1 + z3 - 2 z2
			new complex(z2).zsub(z1).mul(2),           // B = 2 (z2 - z1)
			new complex(z1)                            // C = z1
		);
		this.iff(this.z1, this.z2, abs(this.z1.y) < abs(this.z2.y));
		return this;
	}
	quadinter(z1, z2, z3){// t = this, s = 1-t, result = z1*s*s + 2z2*s*t + z3*t*t
		//return this.spline(z1, z2, z3);
		return this.inter(new complex(this).inter(z1, z2), new complex(this).inter(z2, z3));
	}
	cubiceq(A, B, C, D){// cubic equation solver A z³ + B z² + C z + D = this
	/*
		solve: 
			3 + 5 i = z³ - (2 + 5 i) z² - (5 - 8 i) z + (9 + 2 i)
		code:
			var z = new complex(3, 5).cubiceq(1, [-2, -5], [-5, 8], [9, 2]);
			or
			var z = new complex(3, 5).cubiceq(
				new complex(1),
				new complex(2, 5).neg,
				new complex(5, -8).neg,
				new complex(9, 2)
			);
		roots: 
			z.z1 = 1 + 2 i, z.z2 = 1, z.z3 = 3 i
		check:
			z.asg(z.z1).horner(1, [-2, -5], [-5, 8], [9, 2]); // z = 3 + 5 i
			z.asg(z.z2).horner(1, [-2, -5], [-5, 8], [9, 2]); // z = 3 + 5 i
			z.asg(z.z3).horner(1, [-2, -5], [-5, 8], [9, 2]); // z = 3 + 5 i
	*/
		this.z1 = {}; this.z2 = {}; this.z3 = {}; // result
		var a = new complex(A); // a = A
		if (a.is0){// B z² + C z + D = this
			this.quadraticeq(B, C, D);
			this.z3 = {x: this.z2.x, y: this.z2.y};
		} else {
			var d = new complex(D).zsub(this); // d = D - this
			if (d.is0){// (A z² + B z + C) z = 0
				this.put.zero.obj(this.z3).quadraticeq(A, B, C).pop; // z3 = 0, rest is quadratic, preserve this
			} else {
				const r = new complex(-1, sqrt(3)).div(2); // cis 120°; r³ = 1
				var b = new complex(B); // b = B
				var c = new complex(C); // c = C
				if (b.is0 && c.is0){// A z³ + D = 0
					d.zdiv(a).neg.cbrt
						.obj(this.z1).zmul(r) // z1 = cbrt(-D/A)
						.obj(this.z2).zmul(r) // z2 = z1 rotated 120°
						.obj(this.z3);        // z3 = z2 rotated 120°
				} else {// A z³ + B z² + C z + d = 0
					a.mul(-3); c.zmul(a);               // a = -3 A; c = -3 A C
					var e = new complex(b).sqr;         // e = B²
					var D0 = new complex(e).zadd(c);    // D0 = B² - 3 A C
					e.zmul(b).mul(2); c.zmul(b).mul(3); // e = 2 B³; c = -9 A B C
					var D1 = new complex(a).sqr.zmul(d).mul(3).zadd(e).zadd(c); // D1 = 2 B³ - 9 A B C + 27 A² D
					c.asg(D0).cub.mul(-4).zadd(d.asg(D1).sqr).sqrt.zadd(D1).div(2).cbrt; // c = cbrt((D1 + sqrt(D1² - 4 D0³))/2)
					d.asg(D0).zdiv(c).zadd(c).zadd(b).zdiv(a).obj(this.z1); c.zmul(r); // z1 = (B + c + D0/c)/a; c rotate 120°
					d.asg(D0).zdiv(c).zadd(c).zadd(b).zdiv(a).obj(this.z2); c.zmul(r);
					d.asg(D0).zdiv(c).zadd(c).zadd(b).zdiv(a).obj(this.z3);
				}
			}
		}
		return this;
	}
	beziertimes(z1, z2, z3, z4){// similar to quadratic times but cubic
		// solve (3 (z2 - z3) + z4 - z1) t³ + 3 (z1 + z3 - 2 z2) t² + 3 (z2 - z1) t + z1 = this
		this.cubiceq(
			new complex(z2).zsub(z3).mul(3).zadd(z4).zsub(z1), // A = 3 (z2 - z3) + z4 - z1
			new complex(z2).mul(-2).zadd(z1).zadd(z3).mul(3),  // B = 3 (z1 + z3 - 2 z2)
			new complex(z2).zsub(z1).mul(3),                   // C = 3 (z2 - z1)
			new complex(z1)                                    // D = z1
		);
		this.asg(this.z1)
			.asg(this.z2, abs(this.y) > abs(this.z2.y))
			.asg(this.z3, abs(this.y) > abs(this.z3.y));
		return this;
	}
	bezierinter(z1, z2, z3, z4){// t = this, s = 1-t, result = z1*s*s*s + 3z2*s*s*t + 3z3*s*t*t + z4*t*t*t
		//return this.spline(z1, z2, z3, z4);
		return this.inter(new complex(this).quadinter(z1, z2, z3), new complex(this).quadinter(z2, z3, z4));
	}
	go(z, t = 1){// simplified usage of inter
		return this.asg(new complex(t).inter(this, z));
	}
	halfway(z){
		//return this.go(z, 1/2);
		return this.zadd(z).div(2);
	}
	opposite(z){
		return this.go(z, 2);
	}
	complement(z){
		return this.go(z, 3/2);
	}
	anticomplement(z){
		return this.go(z, 3);
	}
	crossover(z1, z2){
		return this.opposite(new complex(z1).halfway(z2));
	}
	parallel(z1, z2){// this--result is parallel to z1--z2
		return this.zadd(new complex(z2).zsub(z1));
	}
	perp(z1, z2){// this--result is perpendicular to z1--z2
		return this.zadd(new complex(z2).zsub(z1).muli);
	}
	angled(z1, z2, angle){// this--result is slanted by angle to line z1--z2
		return this.asg(new complex(this).ortho(z1, z2).zabout(this, angle));
	}
	intersection(z1, z2, z3, z4){// intersection point of lines z1--z2 and z3--z4
		function cross(p, q){return p.x * q.y - p.y * q.x;}
		var u = new complex(z2).zsub(z1);
		var v = new complex(z4).zsub(z3);
		var d = cross(u, v);
		if (d == 0){// parallel
			return this.inf;
		} else {
			var a = cross(u, z1);
			var b = cross(v, z3);
			u.mul(b); v.mul(a);
			return this.asg(v).zsub(u).div(d);
		}
	}
	bisection(z1, z2, z3, z4 = z1){// intersection of bisectors
		var u = new complex(z1).halfway(z2);
		var v = new complex(z3).halfway(z4);
		return this.intersection(u, new complex(u).perp(z1, z2), v, new complex(v).perp(z3, z4));
	}
	ortho(z1, z2 = {x: 0, y: 0}){// orthogonal projection of this to line z1--z2
		return this.times(z1, z2).real.inter(z1, z2); // 3 : 0 (Jesé 29', Raphaël Varane 56', James Rodríguez 88')
	}
	reflect(z1, z2 = {x: 0, y: 0}){// reflection of this to line z1--z2
		return this.times(z1, z2).conjg.inter(z1, z2);
	}
	dist(x, y){// distance of this to point (x, y)
		return new complex(x, y).zsub(this).abs;
	}
	zdist(z1, z2){// distance to point or line
		if (arguments.length < 2 || new complex(z2).is(z1)) 
			return this.dist(z1.x, z1.y); 
		else
			return this.zdist(new complex(this).ortho(z1, z2));
	}
	azimuth(x, y){// angle of this -- (0, 0) -- (x, y)
		return new complex(x, y).zdiv(this).arg;
	}
	zazimuth(z){
		return azimuth(z.x, z.y);
	}
	linedist(z1, z2 = {x: 0, y: 0}){// directed distance of this to line z1--z2 
		return new complex(z1).zdist(z2) * new complex(this).times(z1, z2).y;
	}
	toward(x, y, len = 1){
		if (this.isEq(x, y))
			return this;
		else
			return this.zadd(new complex(x, y).zsub(this).unit.mul(len));
	}
	ztoward(z, len = 1){// go toward z by given length
		return this.toward(z.x, z.y, len);
	}
	oncircle(circle){// nearest point on circle
		if (this.is(circle))
			return this;
		else
			return this.asg(new complex(circle).ztoward(this, circle.r));
	}
	circlesub(circle){// distance vector from nearest point on circle
		return this.zsub(new complex(this).oncircle(circle));
	}
	circledif(z1, z2, circle){// distance vector from line to circle
		return this.asg(circle).ortho(z1, z2).circlesub(circle);
	}
	isotomicconjg(A, B, C){
		return this.intersection(
			A, new complex().intersection(A, this, B, C).crossover(B, C), 
			B, new complex().intersection(B, this, C, A).crossover(C, A));
	}
	isogonalconjg(A, B, I){
		return this.intersection(
			A, new complex(this).reflect(A, I), 
			B, new complex(this).reflect(B, I));
	}
	harmonicconjg(z1, z2){// z1-------this---z2------result
		// |z1 this| : |this z2| = |z1 result| : |z2 result|
		// result -> (this * (z1 + z2) - 2 * z1 * z2)/(2 * this - (z1 + z2))
		var z = new complex(z1).zadd(z2), w = new complex(2).zmul(z1).zmul(z2);
		// result -> (this * z - w)/(2 * this - z)
		var u = new complex(this).zmul(z).zsub(w), v = new complex(2).zmul(this).zsub(z);
		// result -> u/v
		return this.asg(u).zdiv(v);
	}
	radical(c1, c2){// radical point
		this.asg(c1); var d = this.zdist(c2);
		if (d != 0) this.ztoward(c2, ((sqr(c1.r) - sqr(c2.r))/d + d)/2);
		return this;
	}
	similitude(circle1, circle2){// common tangents intersections
		this.z1 = {x: 1/0, y: 1/0}; 
		this.z2 = {x: 1/0, y: 1/0};
		var u = new complex(circle1).mul(circle2.r);
		var v = new complex(circle2).mul(circle1.r);
		var w = new complex();
		w.asg(v).zadd(u).div(circle1.r + circle2.r).obj(this.z1); // internal
		w.asg(v).zsub(u).div(circle1.r - circle2.r).obj(this.z2); // external
		return this;
	}
	linecircle(point1, point2, circle){// chord, line - circle intersection
		function cross(p, q){return p.x * q.y - p.y * q.x;};
		var u = new complex(point1).zsub(circle);
		var v = new complex(point2).zsub(circle);
		var z = new complex(v).zsub(u);
		var d = z.sqrabs, a = cross(u, v);
		var b = d * circle.r * circle.r - a * a;
		this.z1 = {x: 1/0, y: 1/0}; 
		this.z2 = {x: 1/0, y: 1/0};
		this.put.success = b >= 0;
		if (this.success){
			u.asg(z).mul(a).divi;
			v.asg(z).mul((z.y < 0 ? -1 : 1) * sqrt(b));
			this.asg(u).zadd(v).div(d).zadd(circle).obj(this.z1);
			this.asg(u).zsub(v).div(d).zadd(circle).obj(this.z2);
		}
		return this.pop;
	}
	circlecircle(circle1, circle2){// circle - circle intersection
		if (true){
			this.radical(circle1, circle2);
			return this.linecircle(this, new complex(this).perp(circle1, circle2), circle1);
		} else {
			var z = new complex(circle2).zsub(circle1);
			var d = z.sqrabs, c = sqrt(d);  z.div(c);
			var u = circle1.r, v = circle2.r; u *= u; v *= v;
			var p = d + u - v, q = 4 * d * u - p * p;
			this.z1 = {x: 1/0, y: 1/0}; 
			this.z2 = {x: 1/0, y: 1/0};
			this.put.success = q >= 0;
			if (this.success){
				this.xiy(p, sqrt(q)).div(2 * c);
				this.put.zmul(z).zadd(circle1).obj(this.z1);
				this.pop.conjg.zmul(z).zadd(circle1).obj(this.z2);
			}
			return this.pop;
		}
	}
	circletangent(point, circle){// tangent from point to circle
		if (true){
			var z = new complex(point).halfway(circle); z.r = z.zdist(point);
			this.circlecircle(z, circle);
		} else {
			var z = new complex(point).zsub(circle);
			var d = z.sqrabs, c = sqrt(d);  z.div(c);
			var a = circle.r;
			this.z1 = {x: 1/0, y: 1/0}; 
			this.z2 = {x: 1/0, y: 1/0};
			this.success = c >= a;
			if (this.success){
				var b = sqrt(d - a * a);
				this.xiy(a, b).mul(a).div(c);
				this.put.zmul(z).zadd(circle).obj(this.z1);
				this.pop.conjg.zmul(z).zadd(circle).obj(this.z2);
			}
		}
		this.asg(point); this.r = this.zdist(this.z1);
		return this;
	}
	tangent(circle){// simpler usage
		return this.circletangent(this.z, circle);
	}
	arithSpiral(a, b, t){
		return this.cis(a + b * t, t);
	}
	logSpiral(a, b, t){
		return this.cis(a * exp(b * t), t);
	}
	nautilus(t){
		const a = 1, b = 0.30634896253003312211567570119977; // ln(φ) / (π/2)
		return this.logSpiral(a, b, t);
	}
	seed(n, size = 1, rate = 1, angle = 0){// sunflower seeds
		const f = 2.3999632297286533222315555066336; // π * (3 - sqrt(5))
		return this.zadd(new complex().cis(size * pow(n, rate/2), angle + n * f));
	}
	superellipse(shape = 2, angle, xradius = 1, yradius = xradius, symmetry = 4, u = shape, v = u){
		this.cis(angle * symmetry/4).pos.lcs(xradius, yradius);
		this.xiy(pow(this.x, u), pow(this.y, v));
		return this.cis(pow(this.x + this.y, -1/shape), angle);
	}
	supercircle(shape = 2, angle, radius = 1, symmetry = 4, u = shape, v = u){
		return this.superellipse(shape, angle, radius, radius, symmetry, u, v);
	}
	get cartdev(){// bi-variate uniform distributed complex random
		return this.scl(random(), random());
	}
	get polardev(){// unit circle complex random
		return this.zscl(new complex().cis(sqrt(random()), random(2 * pi)));
	}
	get rectdev(){// unit square complex random
		return this.scl(random(2) - 1, random(2) - 1);
	}
	get expdev(){// bi-variate double sided exponentianl (Laplace) distributed complex random
		function laplace(){return -log(1 - random()) * (random() < 0.5 ? -1 : 1);}
		return this.scl(laplace(), laplace());
	}
	get normaldev(){// bi-variate normal distributied complex random (darts in target)
		return this.zscl(new complex().cis(sqrt(-2 * log(1 - random())), random(2 * pi)));
	}
	get poissondev(){// bi-variate Poisson distributed complex random (result of match)
		function poisson(lambda){
			if (lambda == 0) return 0; else {
				var l = exp(-abs(lambda)), r = -1, p = 1;
				do { r++; p = random(p); } while (p > l);
				if (lambda < 0) r = -r;
				return r;
			};
		}
		return this.xiy(poisson(this.x), poisson(this.y));
	}
	mandelbrot(m = 256){// Mandelbrot set (for testing only)
		// returns 0 to 256; 0: (probably) in set; 256: out of bounds |z|² ≤ 4
		// -2 < x < 1, |y| < 1.2496210676876531737592088948857 for |z|²·|z+1|² ≤ 4
		// eellipse({x: -1/2, y: 0, a: 1.5, b: sqrt((sqrt(17) - 1)/2), o: 0})
		var n = m, z = new complex(this); // |z| ≤ 2 => |z|² ≤ 4
		while (n > 0 && z.sqrabs <= 4) if (z.nop(n--).sqr.zadd(this).is0) n = 0;
		return n;
	}
	trivertex(vertexA, vertexB, vertexC, changed = true){// triangle ABC centers
		var A = new complex(vertexA), B = new complex(vertexB), C = new complex(vertexC);
		this.vertex = {A: {x: A.x, y: A.y}, B: {x: B.x, y: B.y}, C: {x: C.x, y: C.y}};
		this.box = {
			min: {x: min(A.x, B.x, C.x), y: min(A.y, B.y, C.y)}, 
			max: {x: max(A.x, B.x, C.x), y: max(A.y, B.y, C.y)}
		}; 
		this.box.size = {x: this.box.max.x - this.box.min.x, y: this.box.max.y - this.box.min.y};
		this.O = {}; // X3 - circumcircle
		if (changed) {
			this.side = {a: B.zdist(C), b: C.zdist(A), c: A.zdist(B)};
			this.bisection(A, B, C).obj(this.O);
		} else this.obj(this.O);
		this.O.r = this.asg(this.O).zdist(A);
		var a = this.side.a, b = this.side.b, c = this.side.c;
		this.success = (a < b + c) && (b < c + a) && (c < a + b);
		if (this.success){
			var o = A.trap(B) + B.trap(C) + C.trap(A); // signed area
			this.direction = Math.sign(o); // direction of vertices
			var P = a + b + c,  s = P/2,  D = abs(o),  S = 2*D;     // D = Δ
			var abc = a*b*c, aa = a*a, bb = b*b, cc = c*c, ss = s*s;     // common
			var q = aa + bb + cc, p = a*b + b*c + c*a, qq = q*q, pp=p*p; // variables
			var ra = s - a, rb = s - b, rc = s - c; // vertex touch circle radius
			this.vertex.A.r = ra; this.vertex.B.r = rb; this.vertex.C.r = rc;
			this.altitude = {a: S/a, b: S/b, c: S/c}; // altitudes, heights
			this.angle = {A: asin(a/2/this.O.r), B: asin(b/2/this.O.r)}; // angles
			this.angle.C = pi - (this.angle.A + this.angle.B);
			this.perimeter = P; this.area = D; this.semi = s;
			this.omega = atan(4*D, q); // Brocard angle
			this.I = {r: D/s}; // X1 - incircle (weighted average (Aa + Bb + Cc)/(a + b + c))
			this.zero.zadd(new complex(A).mul(a)).zadd(new complex(B).mul(b)).zadd(new complex(C).mul(c)).div(P).obj(this.I);
			this.I.A = {}; this.asg(this.I).ortho(B, C).obj(this.I.A); // incircle
			this.I.B = {}; this.asg(this.I).ortho(C, A).obj(this.I.B); // contact
			this.I.C = {}; this.asg(this.I).ortho(A, B).obj(this.I.C); // points
			this.I.a = this.I.r * sqrt(pp - abc*s - p*ss)/(p - ss); // Adams radius
			this.G = {}; // X2 - centroid 
			this.zero.zadd(A).zadd(B).zadd(C).div(3).obj(this.G);
			this.H = {}; // X4 - orthocenter O----G----+----H
			this.asg(this.O).anticomplement(this.G).obj(this.H); // allways with respect to G
			this.N = {}; // X5 - nine-point circle O-------N-------H
			this.asg(this.O).halfway(this.H).obj(this.N); this.N.r = this.O.r/2;
			this.K = {}; // X6 - symmedian point (isogonal conjugate of G, Lemoine)
			this.asg(this.G).isogonalconjg(A, B, this.I).obj(this.K); // allways with respect to I
			this.Ge = {}; // X7 - Gergonne point (intersection of vertex--contact lines)
			this.intersection(A, this.I.A, B, this.I.B).obj(this.Ge);
			this.Na = {}; // X8 - Nagel point (isotomic conjugate of Ge, anticomplement of I)
			//this.asg(this.Ge).isotomicconjg(A, B, C).obj(this.Na);
			this.asg(this.I).anticomplement(this.G).obj(this.Na);
			this.M = {}; // X9 - Mittentpunkt (complement of Ge) Ge----+----G----M
			this.asg(this.Ge).complement(this.G).obj(this.M);
			this.Sp = {}; // X10 - Spieker center (complement of I)
			this.asg(this.I).complement(this.G).obj(this.Sp); this.Sp.r = this.I.r/2;
			this.F = {}; // X11 - Feuerbach point (common tangent point of incircle and nine-point circle)
			this.asg(this.I).oncircle(this.N).obj(this.F);
			this.X12 = {}; // harmonic conjugate of X11 with respect to X1 and X5
			this.asg(this.F).harmonicconjg(this.I, this.N).obj(this.X12);
			this.L = {}; // X20 - de Longchamps point (Soddy line I--L)
			this.asg(this.H).opposite(this.O).obj(this.L);
			this.shield = {}; // shield circle G----o----H (GH = diameter)
			this.shield.r = this.asg(this.G).halfway(this.H).obj(this.shield).zdist(this.G);
			this.asg(this.K).oncircle(this.shield).zsub(this.G); // ellipse axis inclination vector
			var Z = sqrt((aa*aa + bb*bb + cc*cc) - (aa*bb + bb*cc + cc*aa));// discriminant
			this.Oe = {// Steiner circumellipse
				x: this.G.x, y: this.G.y,
				a: sqrt(q + Z*2)/3, b: sqrt(q - Z*2)/3, c: sqrt(Z) * 2/3, 
				o: this.arg, F1: {}, F2: {}
			}; 
			this.unit.mul(this.Oe.c).zadd(this.Oe).obj(this.Oe.F1).opposite(this.Oe).obj(this.Oe.F2); // foci
			this.Oe.e = this.Oe.c / this.Oe.a;  this.Oe.l = (q - Z*2)/9 / this.Oe.a; // eccentricity, semi-latus rectum
			this.S = {}; // X99 - Steiner point (intersection of circumcircle and circumellipse)
			this.barycentricfun(a, b, c, function(a, b, c){return 1/(b*b - c*c);}).obj(this.S);
			// excircles Ja, Jb, Jc
			this.Ja = {r: D/ra}; this.trilinearxiy(-1,  1,  1).obj(this.Ja);
			this.Jb = {r: D/rb}; this.trilinearxiy( 1, -1,  1).obj(this.Jb);
			this.Jc = {r: D/rc}; this.trilinearxiy( 1,  1, -1).obj(this.Jc);
			// Soddy circles 
			this.SO = {}; // X175 Soddy outer circle (isoperimetric point) exists if a + b + c > 4R + r
			this.SO.r = this.barycentricxiy(a - this.Ja.r, b - this.Jb.r, c - this.Jc.r).obj(this.SO).zdist(A) + ra;
			this.SI = {}; // X176 Soddy inner circle (equal detour point)
			this.SI.r = this.barycentricxiy(a + this.Ja.r, b + this.Jb.r, c + this.Jc.r).obj(this.SI).zdist(A) - ra;
			this.Euler = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Euler line (contains G, H, O, N, L, ...)
			this.Nagel = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Nagel line (contains G, I, Na, Sp, ...)
			this.Soddy = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Soddy line (contains I, Ge, L, SO, SI, ...)
			if (a != b || b != c){ // no lines for equilateral triangle
				var O = this.SO; // circle chords
				this.linecircle(this.G, this.O, O);
				this.asg(this.z1).obj(this.Euler.z1); delete(this.z1);
				this.asg(this.z2).obj(this.Euler.z2); delete(this.z2);
				this.linecircle(this.I, this.G, O);
				this.asg(this.z1).obj(this.Nagel.z1); delete(this.z1);
				this.asg(this.z2).obj(this.Nagel.z2); delete(this.z2);
				this.linecircle(this.I, this.L, O);
				this.asg(this.z1).obj(this.Soddy.z1); delete(this.z1);
				this.asg(this.z2).obj(this.Soddy.z2); delete(this.z2);
			}
		}
		return this.asg(this.O);
	}
	triside(a, b, c, inclination = 0, conjugate = true){// triangle construction from sides
		this.side = {a: a, b: b, c: c};
		this.success = (a < b + c) && (b < c + a) && (c < a + b);
		if (this.success){
			var P = a + b + c,  s = P/2,  D = sqrt(s * (s - a) * (s - b) * (s - c));  
			var A = new complex(0, 0);
			var B = new complex(c, 0);
			var C = new complex(b*b + c*c - a*a, 4*D).div(2*c);
			var Z = new complex().cis(inclination);
			if (conjugate) {C.conjg; Z.conjg;} // vertex C is above or below side c
			var O = new complex().bisection(A, B, C); // circumcenter
			A.zsub(O).zmul(Z).zadd(this);   // translate triangle to (0, 0),
			B.zsub(O).zmul(Z).zadd(this);   // rotate by inclination,
			C.zsub(O).zmul(Z).zadd(this);   // translate to origin,
			this.trivertex(A, B, C, false); // construct from vertices
		}
		return this;
	}
	trieqlat(a, inclination = 0, conjugate = true){// equilateral triangle
		return this.triside(a, a, a, inclination, conjugate);
	}
	trignomon(base, side, inclination = 0, conjugate = true){// gnomon (isosceles)
		return this.triside(side, side, base, inclination, conjugate);
	}
	triright(a, b, inclination = 0, conjugate = true){// right triangle
		return this.triside(a, b, sqrt(a*a + b*b), inclination, conjugate);
	}
	trirightgnomon(a, inclination = 0, conjugate = true){// 45-45-90 triangle
		return this.triright(a, a, inclination, conjugate);
	}
	trimonodrafter(a, inclination = 0, conjugate = true){// 30-60-90 triangle
		return this.triside(a, a * sqrt(3), a * 2, inclination, conjugate);
	}
	trigolden(base, inclination = 0, conjugate = true){// golden triangle
		return this.trignomon(base, base * phi, inclination, conjugate)
	}
	trigoldengnomon(base, inclination = 0, conjugate = true){// golden gnomon
		return this.trignomon(base, base / phi, inclination, conjugate)
	}
	triegypt(base, inclination = 0, conjugate = true){// golden triangle
		return this.triside(base * 3/5, base * 4/5, base, inclination, conjugate)
	}
	triheron(base, inclination = 0, conjugate = true){// Ἥρων ὁ Ἀλεξανδρεύς
		return this.triside(base * 15/14, base * 13/14, base, inclination, conjugate)
	}
	trisrba(base, inclination = 0, conjugate = true){// Србин омиљени троугао
		return this.triside(base * 5/8, base * 7/8, base, inclination, conjugate)
	}
	trikepler(a, inclination = 0, conjugate = true){// Kepler right triangle
		return this.triside(a, a * sqrt(phi), a * phi, inclination, conjugate)
	}
	trikimberling(base, inclination = 0, conjugate = true){// Kimberling golden triangle
		const k = 1.3797865516812012355584834707971; // angles = {φ : φ : 1}, k = 1/(2cos(πφ/(2φ+1))
		return this.trignomon(base, base * k, inclination, conjugate)
	}
	trisilver(a, inclination = 0, conjugate = true){// very acute right triangle
		return this.triright(a, 1, inclination, conjugate)
	}
	trialt(a, b, c, inclination = 0, conjugate = true){// triangle construction from altitudes (heights)
		// harmonic addition: a ÷ b = 1/(1/a + 1/b) = (a * b)/(a + b)
		function o(u, v, w){return (u*v*w)/(u*v + v*w + w*u);} // u ÷ v ÷ w
		this.altitude = {a: a, b: b, c: c};
		var I = o( a,  b,  c); // inradius I
		var A = o(-a,  b,  c); // exradius A
		var B = o( a, -b,  c); // exradius B
		var C = o( a,  b, -c); // exradius C
		this.success = (A > 0) && (B > 0) && (C > 0);
		if (this.success){
			var S = 2 * sqrt(I * A * B * C); // double area
			this.triside(S/a, S/b, S/c, inclination, conjugate); // construct from sides
		}
		return this;
	}
	trilineinc(pointA, pointB, circle){// triangle cnstruction from (tangent) line and incircle (or excircle C)
		var z = new complex().circledif(pointA, pointB, circle); // line-circle distance vector (usually 0)
		var A = new complex(pointA).zsub(z).tangent(circle); // translate line as circle tangent and
		var B = new complex(pointB).zsub(z).tangent(circle); // find tangent points from vertices A and B
		var U = new complex().iff(A.z1, A.z2, // select tangent points
			new complex(A.z1).zdist(A, B) > new complex(A.z2).zdist(A, B));
		var V = new complex().iff(B.z1, B.z2, // which is not lies on side c (A--B)
			new complex(B.z1).zdist(A, B) > new complex(B.z2).zdist(A, B));
		var C = new complex().intersection(A, U, B, V); // find vertex C
		this.success = !C.isInf;
		if (this.success) this.trivertex(A, B, C); // construct from vertices
		return this;
	}
	get trilinear(){// exact trilinear <=> directed distance from sides
		var a = this.direction * this.linedist(this.vertex.B, this.vertex.C);
		var b = this.direction * this.linedist(this.vertex.C, this.vertex.A);
		var c = this.direction * this.linedist(this.vertex.A, this.vertex.B);
		return {a: a, b: b, c: c}
	}
	get trilinearrad(){// normalized in incircle radius
		var t = this.trilinear;
		return {a: t.a/this.I.r, b: t.b/this.I.r, c: t.c/this.I.r}
	}
	get barycentric(){// eaxct barycentric <=> a + b + c = Δ
		var t = this.trilinear;
		var a = t.a * this.side.a/2;
		var b = t.b * this.side.b/2;
		var c = this.area - (a + b);
		return {a: a, b: b, c: c}
	}
	get barycentricnorm(){
		const scale = 1;
		var t = this.barycentric;
		var a = scale * t.a / this.area; 
		var b = scale * t.b / this.area; 
		var c = scale - (a + b);
		return {a: a, b: b, c: c}
	}
	get tripolar(){
		var a = this.zdist(this.vertex.A);
		var b = this.zdist(this.vertex.B);
		var c = this.zdist(this.vertex.C);
		return {a: a, b: b, c: c}
	}
	trilinearxiy(a, b, c){// I = {1 : 1 : 1}
		var k = 2 * this.area / (a * this.side.a + b * this.side.b + c * this.side.c);
		a = a * k - this.I.r; b = b * k - this.I.r; c = c * k - this.I.r;
		var d = cis(this.angle.B);
		this.xiy((a + c * d.x) / d.y, -c * this.direction);
		this.vrotate(this.vertex.B, this.vertex.A).zadd(this.I);
		return this;
	}
	trilinearfun(a, b, c, f){
		return this.trilinearxiy(f(a, b, c), f(b, c, a), f(c, a, b));
	}
	barycentricxiy(a, b, c){// G = {1 : 1 : 1}
		return this.trilinearxiy(a/this.side.a, b/this.side.b, c/this.side.c);
	}
	barycentricfun(a, b, c, f){
		return this.barycentricxiy(f(a, b, c), f(b, c, a), f(c, a, b));
	}
	tripolarxiy(a, b, c){// O = {R : R : R}
		var A = new complex(this.vertex.A); A.r = a;
		var B = new complex(this.vertex.B); B.r = b;
		var C = new complex(this.vertex.C); C.r = c;
		var z = new complex(C).circlecircle(A, B);
		var u = abs(C.r - new complex(z.z1).zdist(C));
		var v = abs(C.r - new complex(z.z2).zdist(C));
		return this.iff(z.z1, z.z2, u < v);
	}
	tripolarfun(a, b, c, f){
		return this.tripolarxiy(f(a, b, c), f(b, c, a), f(c, a, b));
	}
	get stringed(){
		var z = new complex(this).round(100000);
		var x = z.x + '';
		var y = abs(z.y); y = (y == 1 ? '' : y + ' ') + 'і';
		if (z.is0) return '0'; else
		if (z.y == 0) return x; else
		if (z.x == 0) return (z.y < 0 ? '-' : '') + y; else
			return z.x + (z.y < 0 ?' - ':' + ') + y;
	}
	polyvalue(p){
		var z = new complex(this); 
		this.zero;
		for (var i = p.hi; i >= 0; i--)
			this.zmul(z).zadd(new complex(p[i]));
		return this;
	}
	polyprime(p, q){
		q.length = 0;
		for (var i = 1; i < p.length; i++)
			q.push(new complex(p[i]).mul(i).z);
	}
	polydiv(p, q){
		var m = p.hi, n = q.hi, l = p.length - q.length;
		if (l < 0) p = [0]; else {
			var r = new Array(l + 1);
			for (var k = l; k >= 0; k--){
				r[k] = {};
				var z = new complex(p[n + k]).zdiv(new complex(q[n])).obj(r[k]);
				for (var j = n; j >= 0; j--)
					new complex(q[j]).zmul(z).zbus(p[j + k]).obj(p[j + k]);
			}
			p.length = r.length;
			for (var k = 0; k < r.length; k++)
				new complex(r[k]).obj(p[k]);
		}
	}
	polytrim(p){
		while (p.length > 0 && new complex(p[p.hi]).is0) p.length--; // leading zeroes trim
		return p;
	}
	polyarg(p, ...arg){
		var n, i, a;
		p.length = 0; 
		if (arg.length == 1 && Array.isArray(arg[0])) {// only one array parameter
			a = arg[0]; n = a.hi;
			for (i = 0; i <= n; i++){ // incremental copy array as complex polynom
				p.push(new complex(a[i]).z);
				if (p.length == 1) new complex(p[0]).zsub(this).obj(p[0]);
			}
		} else {// zero or more than one parameters or first parameter is not array
			a = arg; n = a.hi;
			for (i = n; i >= 0; i--){ //decremental collect parameters
				p.push(new complex(a[i]).z);
				if (p.length == 1) new complex(p[0]).zsub(this).obj(p[0]);
			}
		}
		return p;
	}
	polysolve(roots, ...arg){// simple Newton polynom roots solver
	/*
		polysolve(roots,  3,  4,  5  )  <=>  3 z² +  4 z    + 5     = this
		polysolve(roots, [3,  4,  5] )  <=>  3    +  4 z    + 5 z²  = this
		polysolve(roots, [3,  4], 5  )  <=> (3    +  4 i) z + 5     = this
		polysolve(roots,  3, [4,  5] )  <=>  3 z  + (4      + 5 i)  = this
		polysolve(roots, [[3, 4], 5] )  <=>  3    +  4 i    + 5 z   = this

		find all polynom roots in complex plane: 
			p[n] zⁿ + p[n - 1] zⁿ⁻¹ + ··· + p[2] z² + p[1] z + p[0] = this
			
		solwe:
			z⁴ - 4 z³ - 19 z² + 46 z + 120 = 0
		code:
			var roots = []; new complex().polysolve(roots, 1, -4, -19, 46, 120);
		roots:
			-3, -2, 4, 5

		solve:
			z⁶ - 6 z⁵ - 26 z⁴ + 144 z³ - 47 z² - 210 z = 0
		code:
			var roots = []; new complex().polysolve(roots, 1, -6, -26, 144, -47, -210, 0);
		roots:
			-5, -1, 0, 2, 3, 7

		solve:
			z⁴ = 1
		code: 
			var roots = []; new complex(1).polysolve(roots, 1, 0, 0, 0, 0);
		roots:
			1, -1, i, -i

		solve:
			z⁴ - (7 + 6 i) z³ - (1 - 30 i) z² + (67 - 4 i) z - (60 + 80 i) = 0
		code:
			var roots = []; new complex().polysolve(roots, 1, [-7, -6], [-1, 30], [67, -4], [-60, -80]);
		roots:
			4, 2 + i, i - 2, 3 + 4 i

		solve:
			z⁴ + z³ + z² + z + 1 = 1
		code:
			var roots = []; new complex(1).polysolve(roots, 1, 1, 1, 1, 1);
		roots:
			0, -1, i, -i
	*/
		roots.length = 0; // reset roots array
		function eps(s, t = {x: 0, y: 0}, e = 1e-17){// precision - Kahan summation
			var a = new complex(s).abs, b = new complex(t).abs, c = a - b, d = (c - a) + b;
			return (abs(c) <= e) && (abs(d) <= e);
		}
		var p = [], n, i; // polynom array, degree and index
		this.polyarg(p, ...arg); this.polytrim(p); // get and trim
		//for (i = p.hi; i >= p.lo; i--) console.log('p[',i,'] =',new complex(p[i]).stringed); // (debug)
		this.del; for (i = p.lo; i <= p.hi; i++) this.put.asg(p[i]).exc; // put polynom into storage (debug)
		i = 0; while (i < p.length && new complex(p[i]).isNum) i++; // check consistency
		
		if (p.length > 1 && i == p.length) {// Solve: p[n]·zⁿ + p[n - 1]·zⁿ⁻¹ + ··· + p[2]·z² + p[1]·z + p[0] = 0
			var u = new complex(),  v = new complex(),  w = new complex();
			var f = new complex(),  g = new complex(),  q = [],  d;
			while (p.length > 1){
				n = p.hi; 
				if (new complex(p[0]).is0) w.zero; else
				switch (n) {
				case 1: // linear
					w.asg(new complex(p[0])).zdiv(new complex(p[1])).neg;
					break;
				case 2: // quadratic
					u.zero.quadraticeq(p[2], p[1], p[0]);
					roots.push(u.z1); roots.push(u.z2);
					p.length = 1; // done
					break;
				case 3: // cubic
					u.zero.cubiceq(p[3], p[2], p[1], p[0]);
					roots.push(u.z1); roots.push(u.z2); roots.push(u.z3);
					p.length = 1; // done
					break;
				default:
					d = 0; for (i = 1; i < n - 1; i++) if (new complex(p[i]).is0) d++; // can be better
					if (d > 0 && d == n - 2) {// p[n] zⁿ + p[0] = 0; z = n-th root of -p[0] / p[n]
						v.asg(new complex(p[0])).zdiv(new complex(p[n])).neg.root(n); // 1st vertex
						for (i = 0; i < n; i++)// roots are vertices of regular central n-polygon
							roots.push(w.asg(v).cyclic(i, n).z);
						p.length = 1; // done
					} else { // Newton part
						i = n * 4 * 1024; // max iterations
						this.polyprime(p, q); // q(z) = p'(z)
						w.xiy(1, 1).polardev; // start with random point in unit circle
						u.inf;  v.inf;
						while (i > 0) {
							i -= 4;
							u.asg(v);  v.asg(w);   // previous 2 iterations
							f.asg(w).polyvalue(p); // evaluate p(w)
							if (f.is0)  break;     // good luck - exact zero
							if (eps(f)) break;     // good enough
							g.asg(w).polyvalue(q); // evaluate p'(w)
							if (g.is0) {           // bad  luck - critical point
								w.asg(1, 1).polardev;
								i += 3; // back 3/4 step
							} else {
								w.zsub(f.zdiv(g)); // w -= p(w) / p'(w)
								if (eps(w, u) || eps(w, v)) // almost done
								if (i > 20) i = 20; // few more iterations
							}
						}
					}
				}
				if (p.length > 1){
					roots.push(w.z); // put root in list
					this.polydiv(p, [w.neg.z, 1]); // p /= (z - w) removes root w from polynom p
				}
			}
		}

		return this;
	}
	get sto(){// storage status
		return typeof this.mem === 'undefined' ? -1 : this.mem.length;
	}
	get del(){// delete storage
		if (this.sto >= 0) delete(this.mem);
		return this;
	}
	get clr(){// clear last from storage
		if (this.sto > 0) {
			this.mem.pop();
			if (this.mem.length == 0) this.del;
		}
		return this;
	}
	get put(){// put into storage
		if (this.sto < 0) this.mem = [];
		this.mem.push(this.z);
		return this;
	}
	get pop(){// get and clear
		return this.get.clr;
	}
	get set(){// set storage
		if (this.sto > 0) this.obj(this.mem[this.mem.hi]);
		return this;
	}
	get get(){// get from storage
		if (this.sto > 0) this.asg(this.mem[this.mem.hi]);
		return this;
	}
	get exc(){// exchange with storage
		if (this.sto > 0) this.swap(this.mem[this.mem.hi]);
		return this;
	}
	get rot(){// rotate storage
		if (this.sto > 0) this.mem.push(this.mem.shift());
		return this;
	}
	get rol(){// rotate and get
		return this.rot.get;
	}
}

// Не могу више да куцам Math. Ја сам паскал програмер.
min = Math.min, max = Math.max;
abs = function(z){return z instanceof Object ? Math.hypot(z.x, z.y) : Math.abs(z)};
sqrt = Math.sqrt, sqr = function(x){return x * x;};
exp = Math.exp, log = Math.log; pow = Math.pow;
sin = Math.sin, cos = Math.cos, tan = Math.tan;
asin = Math.asin, acos = Math.acos;
atan = function(s, c = 1){return s instanceof Object ? Math.atan2(s.y, s.x) : Math.atan2(s, c);}
cis = function(a, x = 1, y = x){return {x: x * cos(a), y: y * sin(a)};}
pi = Math.PI;          // π = 3.1415926535897932384626433832795
phi = (sqrt(5) + 1)/2; // φ = 1.6180339887498948482045868343656
random = function(x = 1){return x * Math.random();}


Object.defineProperty(Array.prototype, 'lo', {
	enumerable: false, configurable: false,
	get() { return 0; }
});

Object.defineProperty(Array.prototype, 'hi', {
	enumerable: false, configurable: false,
	get() { return this.length - 1; }
});


Object.defineProperty(CanvasRenderingContext2D.prototype, 'begin', {
	get() { this.beginPath(); return this; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'close', {
	get() { this.closePath(); return this; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'width', {
	enumerable: false, configurable: false,
	get() { return this.canvas.width; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'height', {
	enumerable: false, configurable: false,
	get() { return this.canvas.height; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'sizeX', {
	enumerable: false, configurable: false,
	get() { return this.width; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'sizeY', {
	enumerable: false, configurable: false,
	get() { return this.height; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'size', {
	enumerable: false, configurable: false,
	get() { return {x: this.sizeX, y: this.sizeY}; }
});

CanvasRenderingContext2D.prototype.clear = function(){
	this.save();
	this.clearRect(0, 0, this.sizeX, this.sizeY);
	this.restore();
	return this;
}

CanvasRenderingContext2D.prototype.strokefill = function () {
	this.stroke(); this.fill(); return this;
}

CanvasRenderingContext2D.prototype.zmoveTo = function(z){
	this.moveTo(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zlineTo = function(z){
	this.lineTo(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zline = function(z1, z2){
	this.zmoveTo(z1).zlineTo(z2); return this;
}

CanvasRenderingContext2D.prototype.zrect = function(z1, z2){
	var z = new complex(z2).zsub(z1);
	this.rect(z1.x, z1.y, z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zfillRect = function(z1, z2){
	var z = new complex(z2).zsub(z1);
	this.fillRect(z1.x, z1.y, z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zstrokeRect = function(z1, z2){
	var z = new complex(z2).zsub(z1);
	this.strokeRect(z1.x, z1.y, z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zclearRect = function(z1, z2){
	var z = new complex(z2).zsub(z1);
	this.clearRect(z1.x, z1.y, z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zquadraticCurveTo = function(z1, z2){
	this.quadraticCurveTo(z1.x, z1.y, z2.x, z2.y); return this;
}

CanvasRenderingContext2D.prototype.zbezierCurveTo = function(z1, z2, z3){
	this.bezierCurveTo(z1.x, z1.y, z2.x, z2.y, z3.x, z3.y); return this;
}

CanvasRenderingContext2D.prototype.zarc = function(z, r, s, e, c = false){
	if (r < 0){r = -r; c = !c;}
	this.arc(z.x, z.y, r, s, e, c); return this;
}

CanvasRenderingContext2D.prototype.carc = function(z, s, e, c = false){
	return this.zarc(z, z.r, s, e, c);
}

CanvasRenderingContext2D.prototype.zarcTo = function(z1, z2, r){
	this.arcTo(z1.x, z1.y, z2.x, z2.y, Math.abs(r)); return this;
}

CanvasRenderingContext2D.prototype.isZPointInPath = function(z){// must be improved
	return this.isPointInPath(z.x, z.y);
}

CanvasRenderingContext2D.prototype.isZPointInStroke = function(z){// must be improved
	return this.isPointInStroke(z.x, z.y);
}

CanvasRenderingContext2D.prototype.zscale = function(z){
	this.scale(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zrotate = function(z){
	this.rotate(Math.atan2(z.y, z.x)); return this;
}

CanvasRenderingContext2D.prototype.ztranslate = function(z){
	this.translate(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zsetTransform = function(scale, skew, move){
	this.setTransform(scale.x, skew.x, skew.y, scale.y, move.x, move.y); return this;
}


CanvasRenderingContext2D.prototype.circle = function (x, y, r) {
	this.arc(x, y, Math.abs(r), 0, 2 * Math.PI); return this;
}

CanvasRenderingContext2D.prototype.zcircle = function (z, r) {
	return this.circle(z.x, z.y, r);
}

CanvasRenderingContext2D.prototype.ccircle = function (c) {
	return this.zcircle(c, c.r);
}

CanvasRenderingContext2D.prototype.zellipse = 
function(z, a, b, rot = 0, s = 0, e = 2*Math.PI, c = false){
	this.ellipse(z.x, z.y, Math.abs(a), Math.abs(b), rot, s, e, c); return this;
}

CanvasRenderingContext2D.prototype.eellipse = 
function(z, s = 0, e = 2*Math.PI, c = false){
	return this.zellipse(z, z.a, z.b, z.o, s, e, c);
}

CanvasRenderingContext2D.prototype.fellipse = 
function(F1x, F1y, F2x, F2y, b, s = 0, e = 2*Math.PI, c = false){
	var dx = F1x - F2x, dy = F1y - F2y, x = F1x - dx/2, y = F1y - dy/2; o = Math.atan2(dy, dx);
	var d = dx*dx + dy*dy, a = Math.sqrt(b*b + d); b = Math.abs(b);
	this.ellipse(x, y, a, b, o, s, e, c); return this;
}


CanvasRenderingContext2D.prototype.superellipse = 
function(shape = 2, xpos, ypos, xradius = 1, yradius = xradius, azimuth = 0, symmetry = 4, u = shape, v = u){
	var steps = symmetry * Math.max(9, xradius + yradius);
	if (steps > 0){
		var z = new complex(), theta;
		for (var i = 0; i < steps; i++){
			theta = 2 * i * Math.PI/steps;
			z.superellipse(shape, theta, xradius, yradius, symmetry, u, v).rotate(azimuth).add(xpos, ypos);
			if (i == 0)
				this.zmoveTo(z);
			else
				this.zlineTo(z); // to do: can be improved with spline and less steps
		}
		this.close;
	}
	return this;
}

CanvasRenderingContext2D.prototype.zsuperellipse = 
function(shape = 2, z, xradius = 1, yradius = xradius, azimuth = 0, symmetry = 4, u = shape, v = u){
	return this.superellipse(shape, z.x, z.y, xradius, yradius, azimuth, symmetry, u, v);
}

CanvasRenderingContext2D.prototype.esuperellipse = 
function(shape = 2, z, symmetry = 4, u = shape, v = u){
	return this.zsuperellipse(shape, z, z.a, z.b, z.o, symmetry, u, v);
}

CanvasRenderingContext2D.prototype.supercircle = 
function(shape = 2, xpos, ypos, radius = 1, azimuth = 0, symmetry = 4, u = shape, v = u){
	return this.superellipse(shape, xpos, ypos, radius, radius, azimuth, symmetry, u, v);
}

CanvasRenderingContext2D.prototype.zsupercircle = 
function(shape = 2, z, radius = 1, azimuth = 0){
	return this.supercircle(shape, z.x, z.y, radius, azimuth);
}

CanvasRenderingContext2D.prototype.csupercircle = 
function(shape = 2, z, azimuth = 0){
	return this.zsupercircle(shape, z, z.r, azimuth);
}

CanvasRenderingContext2D.prototype.supertable = 
function(shape = 2, xpos, ypos, radius = 1, azimuth = 0, symmetry = 4, u = shape, v = u){
	return this.superellipse(shape, xpos, ypos, radius, radius * (Math.sqrt(5) - 1)/2, azimuth, symmetry, u, v);
}

CanvasRenderingContext2D.prototype.triangle = function(t) {
	if (t.success){
		return this
		.zmoveTo(t.vertex.A)
		.zlineTo(t.vertex.B)
		.zlineTo(t.vertex.C)
		.close;
	} else return this;
}

CanvasRenderingContext2D.prototype.roundRect = function(x, y, w, h, r = 0) {
	// https://stackoverflow.com/questions/1255512/
	// https://stackoverflow.com/users/167531/grumdrig
	if (w < 2 * r) r = w / 2;
	if (h < 2 * r) r = h / 2;
	this.moveTo(x+r, y);
	this.arcTo(x+w, y,   x+w, y+h, r);
	this.arcTo(x+w, y+h, x,   y+h, r);
	this.arcTo(x,   y+h, x,   y,   r);
	this.arcTo(x,   y,   x+w, y,   r);
	return this.close;
}

CanvasRenderingContext2D.prototype.zroundrect = function (z1, z2, r = 0) {
	var z = new complex(z2).zsub(z1);
	return this.roundRect(z1.x, z1.y, z.x, z.y, r);
}

CanvasRenderingContext2D.prototype.rgba = function(r, g, b, a = 1){
	return 'rgba('+r+', '+g+', '+b+', '+a+') ';
}

CanvasRenderingContext2D.prototype.shadowOffsetZ = function(z){
	this.shadowOffsetX = z.x;
	this.shadowOffsetY = z.y;
	return this;
}

CanvasRenderingContext2D.prototype.zfillText = function(text, z, width = -1){
	if (width < 0)
		this.fillText(text, z.x, z.y);
	else
		this.fillText(text, z.x, z.y, w);
	return this;
}

CanvasRenderingContext2D.prototype.zstrokeText = function(text, z, width = -1){
	if (width < 0)
		this.strokeText(text, z.x, z.y);
	else
		this.strokeText(text, z.x, z.y, w);
	return this;
}

// Проба
CanvasRenderingContext2D.prototype.TextL = function(z, text){
	var m = this.measureText(text);
	this.zfillText(text, z);
	this.last = {l: {x: z.x, y: z.y}, r: {x: z.x + m.width, y: z.y}, w: m.width};
	return this;
}

CanvasRenderingContext2D.prototype.TextR = function(z, text){
	var m = this.measureText(text);
	return this.TextL(new complex(z).sub(m.width), text);
}

CanvasRenderingContext2D.prototype.TextC = function(z, text){
	var m = this.measureText(text);
	return this.TextL(new complex(z).sub(m.width/2), text);
}
