class complex {
/*
	Complex arithmetic with canvas extension
	
	JavaScript implementation by Србислав Д. Нешић, November 2019
	srbislav.nesic@gmail.com
	Nasumica Agencija
	
	https://github.com/Nasumica/JScomplex

	Euler identity e ^ (i * π) + 1 = 0: new complex(pi).muli.exp.add(1);
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
		if (u instanceof Object) {
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
		var r = pi/2; // few exact values of right angle
		if (theta == 0  || theta == tau)    this.xiy(+rho, 0); else
		if (theta == r  || theta == -3 * r) this.xiy(0, +rho); else 
		if (theta == -r || theta == 3 * r)  this.xiy(0, -rho); else 
		if (theta == pi || theta == -pi)    this.xiy(-rho, 0); else
			this.xiy(rho * cos(theta), rho * sin(theta));
		return this;
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
	get cab(){// taxicab distance
		return abs(this.x) + abs(this.y);
	}
	get sqrabs(){// |z|² = z * z' = ρ²
		return this.x * this.x + this.y * this.y;
	}
	get abs(){// magnitude ρ
		return this.x == 0 || this.y == 0 ? this.cab : sqrt(this.sqrabs);
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
	get dbl(){
		return this.scl(2);
	}
	lcs(x, y = x){
		return this.xiy(this.x / x, this.y / y);
	}
	zlcs(z){// scale down
		return this.lsc(z.x, z.y);
	}
	get half(){
		return this.lcs(2);
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
		if (this.isInf) return this.zero; else
		if (this.y == 0) return this.xiy(1 / this.x); else 
		if (this.x == 0) return this.xiy(1 / this.y).divi; else
			return this.conjg.lcs(this.sqrabs);
	}
	mod(x = 1, y = 0){
		var z = new complex(this).div(x, y).floor().mul(x, y);
		return this.zsub(z);
	}
	zmod(z){
		return this.mod(z.x, z.y);
	}
	hadd(x, y = 0){
		if (this.isInf) this.xiy(x, y); else {
			var z = new complex(x, y);
			if (z.is0) this.zero; else
			if (!z.isInf) {
				z.zadd(this);
				this.mul(x, y).zdiv(z);
			}
		}
		return this;
	}
	zhadd(z){// harmonic addition
		return this.hadd(z.x, z.y);
	}
	hsub(x, y = 0){
		return this.hadd(-x, -y);
	}
	zhsub(z){// harmonic subtraction
		return this.hsub(z.x, z.y);
	}
	hmul(x, y = 0){
		return this.div(x, y);
	}
	zhmul(z){// harmonc multiplication
		return this.hmul(z.x, z.y);
	}
	hdiv(x, y = 0){
		return this.mul(x, y);
	}
	zhdiv(z){// harmonic division
		return this.hdiv(z.x, z.y);
	}
	get sqr(){// this²
		return this.xiy((this.x + this.y) * (this.x - this.y), 2 * this.x * this.y);
		//return this.zmul(this);
	}
	get cub(){// this³
		var xx = this.x * this.x, yy = this.y * this.y;
		return this.xiy(this.x * (xx - 3*yy), this.y * (3*xx - yy));
		//return this.zmul(new complex(this).sqr);
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
				return this.xiy(Math.sign(this.x) * pow(abs(this.x), 1/3)); else
			if (this.x == 0)
				return this.xiy(Math.sign(this.y) * pow(abs(this.y), 1/3)).divi; else
			return this.root(3);
	}
	get unit(){// unit vector, sign
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
	get inc(){
		this.x++; return this;
	}
	get dec(){
		this.x--; return this;
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
	get compmod(){// complementary modulus = sqrt(1 - z²)
		return this.sqr.neg.inc.sqrt;
	}
	pagm(i = 32){// pairwise arithmetic-geometric mean
		var z = new complex();
		while (i-- > 0 && this.x != this.y && !this.is(z))
			this.obj(z).xiy((this.x + this.y)/2, sqrt(this.x * this.y));
		return this;
	}
	agm(z, i = 32){// arithmetic-geometric mean
		var u = new complex(z), v = new complex();
		while (i-- > 0 && !this.is(u) && !this.is(v))
			this.obj(v).halfway(u).nop(u.zmul(v).sqrt);
		return this;
	}
	get ellipticK(){// complete elliptic integral of the first kind, K(0) = π/2
		// Mathematica function EllipticK[z] = z.sqrt.ellipticK
		var z = {}; return this.inc.obj(z).vid(2).dec.agm(1).dbl.zmul(z).vid(pi);
	}
	pap(g = 9.80665){// pendulum amplitude period: ρ = length, θ = angle
		return 4 * sqrt(this.abs/g) * new complex(sin(this.arg/2)).ellipticK.x;
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
		if (y == 0 && x ==  2) return this.sqr; else
		if (y == 0 && x ==  3) return this.cub; else
		if (y == 0 && x > 0 && this.is0) return this; else
			return this.log.mul(x, y).exp;
	}
	zpow(z){
		return this.pow(z.x, z.y);
	}
	npow(n){// thisⁿ (integer n); De Moivre's formula: (cis θ)ⁿ = cis (θ·n)
		var z = new complex(this), i = Math.floor(abs(n));
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
		if (y == 0 && x ==  2) return this.sqrt; else
		if (y == 0 && x > 0 && this.is0) return this; else
			return this.log.div(x, y).exp;
	}
	zroot(z){
		return this.root(z.x, z.y);
	}
	get pow2(){// 2^z (must improve)
		return this.lcs(1.442695040888963407359924681001892137427).exp;
	}
	get powpi(){// pi^z (must improve)
		return this.scl(1.144729885849400174143427351353058711647).exp;
	}
	get pow10(){// 10^z (must improve)
		return this.scl(2.302585092994045684017991454684364207601).exp;
	}
	get sinh(){// (e^z - e^-z)/2
		return this.exp.zsub(new complex(this).recip).half;
	}
	get cosh(){// (e^z + e^-z)/2
		return this.exp.zadd(new complex(this).recip).half;
	}
	get tanh(){// (e^2z - 1)/(e^2z + 1)
		this.dbl.exp; var z = new complex(this);
		return this.dec.zdiv(z.inc);
	}
	get csch(){// cosecant
		return this.sinh.recip;
	}
	get sech(){// secant
		return this.cosh.recip;
	}
	get coth(){// cotangent
		return this.tanh.recip;
	}
	get asinh(){// ln(z + sqrt(z² + 1))
		return this.zadd(new complex(this).sqr.inc.sqrt).log;
	}
	get acosh(){// ln(z + sqrt(z² - 1))
		return this.zadd(new complex(this).sqr.dec.sqrt).log;
	}
	get atanh(){// ln((1 + z)/(1 - z))/2
	//  (1 + z)/(1 - z) = 
	//  (z + 1)/(1 - z) =
	// -(z + 1)/(z - 1) = 
	// -(2 + z - 1)/(z - 1) = 
	// -2/(z - 1) - (z - 1)/(z - 1) = 
	// -2/(z - 1) - 1
		return this.dec.recip.dbl.neg.dec.log.half; // :)
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
	get csc(){// cosecant
		return this.sin.recip;
	}
	get sec(){// secant
		return this.cos.recip;
	}
	get cot(){// cotangent
		return this.tan.recip;
	}
	get hav(){// haversine
		return this.half.sin.sqr;
	}
	get crd(){// chord
		return this.half.sin.dbl;
	}
	get sinc(){// sin(z)/z
		var z = {}; return this.is0 ? this.one : this.obj(z).sin.zdiv(z);
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
	get ahav(){
		return this.sqrt.asin.dbl;
	}
	get acrd(){
		return this.half.asin.dbl;
	}
	geodist(x, y, r = 6371008.8){// mean Earth radius r = 6371008.8 m
		return r * ahav(hav(x - this.x) + hav(y - this.y) * cos(x) * cos(this.x));
	}
	zgeodist(z, r = 6371008.8){
		return geodist(z.x, z.y, r);
	}
	horner(){// Horner's scheme polynom evaluate
		var z = {}; this.obj(z).zero; // z = this; this = 0;
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
		return this.zmul(new complex().cis(angle));
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
		return this.rotate(tau * m/n);
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
			var d = new complex(b.neg).sqr.zsub(c.dbl.zmul(a.dbl)).sqrt;
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
				d.obj(this.z3); // z3 = 0
				this.quadraticeq(A, B, C); // rest is quadratic
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
		this.iff(this.z1, this.z2, abs(this.z1.y) < abs(this.z2.y))
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
	halfway(z){// ○-----X-----Z
		//return this.go(z, 1/2);
		return this.zadd(z).half;
	}
	opposite(z){// ○-----Z-----X
		return this.go(z, 2);
	}
	complement(z){// ○-----+-----Z-----X
		return this.go(z, 3/2);
	}
	anticomplement(z){// X-----+-----Z-----o
		return this.go(z, 3);
	}
	get metalic(){// metalic mean = this + 1/(this + 1/(this + 1/(this + ...)))
		return this.halfway(new complex(this).sqr.add(4).sqrt); // x = y - 1/y
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
	angled(z1, z2, angle = 0){// direction to line by given angle
		return this.asg(new complex(this).ortho(z1, z2).zabout(this, angle));
	}
	intersection(z1, z2, z3, z4){// line (z1--z2) - line (z3--z4) intersection point
		function cross(p, q){return p.x * q.y - p.y * q.x;} // p × q
		var u = new complex(z2).zsub(z1); // line1 distance vector z2 - z1
		var v = new complex(z4).zsub(z3); // line2 distance vector z4 - z3
		var d = cross(u, v); // vectors cross product d = u × v
		if (d == 0){// parallel
			this.inf;
		} else {                                      // (u × z1) * v - (v × z3) * u
			var a = cross(u, z1), b = cross(v, z3);   // ---------------------------
			this.asg(v.mul(a)).zsub(u.mul(b)).div(d); //           u × v
		}
		return this;
	}
	bisection(z1, z2, z3, z4 = z1){// intersection of bisectors
		var u = new complex(z1).halfway(z2);
		var v = new complex(z3).halfway(z4);
		return this.intersection(u.z, u.perp(z1, z2), v.z, v.perp(z3, z4));
	}
	ortho(z1, z2){// orthogonal projection of this to line z1--z2
		if (z1.x == z2.x && z1.y == z2.y)
			this.asg(z1);
		else
			this.times(z1, z2).real.inter(z1, z2); // 3 : 0 (Jesé 29', Raphaël Varane 56', James Rodríguez 88')
		return this;
	}
	reflect(z1, z2){// reflection of this to line z1--z2
		if (z1.x == z2.x && z1.y == z2.y)
			this.opposite(z1);
		else
			this.times(z1, z2).conjg.inter(z1, z2);
		return this;
	}
	dist(x, y = 0){// absolute distance of this to point (x, y)
		return new complex(x, y).zsub(this).abs;
	}
	zdist(z1, z2){// absolute distance to point or line
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
	linedist(z1, z2){// directed distance of this to line z1--z2 
		return new complex(z1).zdist(z2) * new complex(this).times(z1, z2).y;
	}
	toward(x, y, len = 1){
		if (!this.isEq(x, y)) this.zadd(new complex(x, y).zsub(this).unit.mul(len));
		return this;
	}
	ztoward(z, len = 1){// go toward z by given length
		return this.toward(z.x, z.y, len);
	}
	slide(z1, z2, len = 1){// slide parallel to z1--z2
		return this.zadd(new complex(z2).zsub(z1).unit.mul(len));
	}
	oncircle(circle){// nearest point on circle
		if (!this.is(circle)) this.asg(new complex(circle).ztoward(this, circle.r));
		return this;
	}
	circledist(circle){// distance vector from nearest point on circle
		return this.zsub(new complex(this).oncircle(circle));
	}
	circledif(z1, z2, circle){// distance vector from line to circle
		return this.asg(circle).ortho(z1, z2).circledist(circle);
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
	radical(circle1, circle2){// radical point
		this.asg(circle1); var d = this.zdist(circle2); // centers distance
		if (d != 0) {
			d = ((sqr(circle1.r) - sqr(circle2.r))/d + d)/2; // distance from C1 center
			this.ztoward(circle2, d); // run
		}
		return this;
	}
	similitude(circle1, circle2){// common tangents intersections
		this.z1 = {x: 1/0, y: 1/0}; this.z2 = {x: 1/0, y: 1/0}; // homothetic centers
		var u = new complex(circle1).mul(circle2.r);
		var v = new complex(circle2).mul(circle1.r);
		new complex()
			.asg(v).zadd(u).div(circle1.r + circle2.r).obj(this.z1)  // internal
			.asg(v).zsub(u).div(circle1.r - circle2.r).obj(this.z2); // external
		return this;
	}
	siminter(circle1, circle2){// internal similitude point
		return this.asg(new complex().similitude(circle1, circle2).z1);
	}
	simouter(circle1, circle2){// external similitude point
		return this.asg(new complex().similitude(circle1, circle2).z2);
	}
	chord(point1, point2, circle){// line - circle intersection
	/*
		                      Z₁
		                r   / |
		                  /   | sqrt(r² - d²)
		                /     |
		(-------------○-------Z----)
		                  d   |
		                      |
		                      |
		                      Z₂
	*/
		this.z1 = {x: 1/0, y: 1/0}; this.z2 = {x: 1/0, y: 1/0};
		var z = new complex(circle).ortho(point1, point2);          // project center to line
		var c = sqr(circle.r) - new complex(z).zsub(circle).sqrabs; // |OZ₁| = r, |OZ| = d, c = r² - d² = |ZZ₁|²
		if (c > 0)                                                  // Z must be in circle else not intersect
			new complex(point1).zsub(point2).unit.mul(sqrt(c))      // ZZ₁
				.zadd(z).obj(this.z1)                               // Z₁
				.opposite(z).obj(this.z2);                          // Z₂
		return this;
	}
	lens(circle1, circle2){// common chord, circle - circle intersection
	/*
		                     |
		(----------○₁----(---R--)------○₂-------------)
		                     |
	*/
		// points of concurrence lies on radical axis
		var r = new complex().radical(circle1, circle2); 
		return this.chord(r.z, r.perp(circle1, circle2), circle1);
	}
	tangent(circle){// tangent from point to circle
	/*
		(--------(○--------)---C------------P)
	*/
		// C is circle with diameter OP; ∟OZ₁P = ∟OZ₂P = 90°
		var c = new complex(this).halfway(circle); c.r = c.zdist(circle);
		return this.lens(c, circle);
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
		return this
			.cis(angle * symmetry/4).pos.lcs(xradius, yradius)
			.xiy(pow(this.x, u), pow(this.y, v))
			.cis(pow(this.x + this.y, -1/shape), angle);
	}
	supercircle(shape = 2, angle, radius = 1, symmetry = 4, u = shape, v = u){
		return this.superellipse(shape, angle, radius, radius, symmetry, u, v);
	}
	get cartdev(){// uniform distributed complex random
		return this.scl(random(), random());
	}
	get polardev(){// unit circle complex random
		return this.zscl(cis(random(tau), sqrt(random())));
	}
	get rectdev(){// unit square complex random
		return this.scl(random(-1, 1), random(-1, 1));
	}
	get expdev(){// bi-variate double sided exponentianl (Laplace) distributed random
		function laplace(){return -log(1 - random()) * (random() < 0.5 ? -1 : 1);}
		return this.scl(laplace(), laplace());
	}
	get normaldev(){// bi-variate normal distributied random (darts in target)
		return this.zscl(cis(random(tau), sqrt(-2 * log(1 - random()))));
	}
	get poissondev(){// bi-variate Poisson distributed random (result of match)
		function poisson(lambda){
			if (lambda == 0) return 0; else {
				var l = exp(-abs(lambda)), r = -1, p = 1;
				do { r++; p *= random(); } while (p > l);
				if (lambda < 0) r = -r;
				return r;
			};
		}
		return this.xiy(poisson(this.x), poisson(this.y));
	}
	mandelbrot(m = 1000, r = 2){// Mandelbrot set (m = max number of iterations)
		var z = new complex(this), w = {}; // |z| ≤ 2 => |z|² ≤ 4 (avoid sqrt)
		while (m > 1 && z.obj(w).sqrabs <= 4) if (z.pow(r).zadd(this).nop(m--).is(w)) m = 0;
		// 0: in set; 1: probably in set; 2 to 999: not in set; 1000: out of bounds
		return m;
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
			var P = a + b + c,  s = P/2,  D = abs(o),  S = 2*D; // D = Δ
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
			this.barycentricfun(aa, bb, cc, function(a, b, c){return 1/(b - c);}).obj(this.S);
			// excircles Ja, Jb, Jc
			this.Ja = {r: D/ra}; this.trilinearxiy(-1,  1,  1).obj(this.Ja);
			this.Jb = {r: D/rb}; this.trilinearxiy( 1, -1,  1).obj(this.Jb);
			this.Jc = {r: D/rc}; this.trilinearxiy( 1,  1, -1).obj(this.Jc);
			// Soddy circles 
			this.SO = {}; // X175 Soddy outer circle (isoperimetric point) exists if a + b + c > 4R + r
			this.SO.r = this.barycentricxiy(a - this.Ja.r, b - this.Jb.r, c - this.Jc.r).obj(this.SO).zdist(A) + ra;
			this.SI = {}; // X176 Soddy inner circle (equal detour point)
			this.SI.r = this.barycentricxiy(a + this.Ja.r, b + this.Jb.r, c + this.Jc.r).obj(this.SI).zdist(A) - ra;
			this.Euler = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Euler line (contains O, G, H, N, L, ...)
			this.Nagel = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Nagel line (contains G, I, Na, Sp, ...)
			this.Soddy = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Soddy line (contains I, L, Ge, SO, SI, ...)
			if (a != b || b != c){ // no lines for equilateral triangle
				var O = this.SO; // lines as this circle chords
				this.chord(this.O, this.G, O) // O--G
					.asg(this.z1).obj(this.Euler.z1)
					.asg(this.z2).obj(this.Euler.z2);
				this.chord(this.G, this.I, O) // G--I
					.asg(this.z1).obj(this.Nagel.z1)
					.asg(this.z2).obj(this.Nagel.z2);
				this.chord(this.I, this.L, O) // I--L
					.asg(this.z1).obj(this.Soddy.z1)
					.asg(this.z2).obj(this.Soddy.z2);
				delete(this.z1); delete(this.z2);
			}
		}
		return this.asg(this.O);
	}
	trisss(a, b, c, inclination = 0, conjugate = true){// side side side
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
	trisas(c, A, b, inclination = 0, conjugate = true){// side angle side
		return this.trisss(new complex().cis(b, A).dist(c), b, c, inclination, conjugate); 
	}
	trissa(c, b, U, inclination = 0, conjugate = true){// side side angle
		var u = max(b, c), v = min(b, c);
		return this.trisas(u, pi - (asin(sin(U) * v/u) + U), v, inclination, conjugate); 
	}
	triasa(A, c, B, inclination = 0, conjugate = true){// angle side angle
		var C = new complex().intersection(
			{x: 0, y: 0}, new complex(c, 0).about(0, 0, +A),
			{x: c, y: 0}, new complex(0, 0).about(c, 0, -B));
		return this.trisss(C.dist(c, 0), C.dist(0, 0), c, inclination, conjugate);
	}
	trieqlat(a, inclination = 0, conjugate = true){// equilateral triangle
		return this.trisss(a, a, a, inclination, conjugate);
	}
	triiso(base, side, inclination = 0, conjugate = true){// isosceles
		return this.trisss(side, side, base, inclination, conjugate);
	}
	trignomon(k, base, inclination = 0, conjugate = true){// gnomon
		return this.triiso(base, base * k, inclination, conjugate);
	}
	triright(a, b, inclination = 0, conjugate = true){// right triangle
		return this.trisss(a, b, sqrt(a*a + b*b), inclination, conjugate);
	}
	trirightgnomon(a, inclination = 0, conjugate = true){// 45-45-90 triangle
		return this.triright(a, a, inclination, conjugate);
	}
	trimonodrafter(a, inclination = 0, conjugate = true){// 30-60-90 triangle
		return this.trisss(a, a * sqrt(3), a * 2, inclination, conjugate);
	}
	trimetalic(m, base, inclination = 0, conjugate = true){// metalic triangle
		return this.trignomon((m + Math.hypot(m, 2))/2, base, inclination, conjugate)
	}
	trigolden(base, inclination = 0, conjugate = true){// golden triangle
		return this.trimetalic(1, base, inclination, conjugate)
	}
	trisilver(base, inclination = 0, conjugate = true){// silver triangle
		return this.trimetalic(2, base, inclination, conjugate)
	}
	tribronze(base, inclination = 0, conjugate = true){// bronze triangle
		return this.trimetalic(3, base, inclination, conjugate)
	}
	trigoldengnomon(base, inclination = 0, conjugate = true){// golden gnomon
		return this.trimetalic(-1, base, inclination, conjugate)
	}
	triegypt(base, inclination = 0, conjugate = true){// 3-4-5 triangle
		return this.trisss(base * 3/5, base * 4/5, base, inclination, conjugate)
	}
	triheron(base, inclination = 0, conjugate = true){// Ἥρων ὁ Ἀλεξανδρεύς
		return this.trisss(base * 15/14, base * 13/14, base, inclination, conjugate)
	}
	trisrba(base, inclination = 0, conjugate = true){// мој омиљени троугао
		return this.trisss(base * 5/8, base * 7/8, base, inclination, conjugate)
	}
	trikepler(a, inclination = 0, conjugate = true){// Kepler right triangle
		return this.trisss(a, a * sqrt(phi), a * phi, inclination, conjugate)
	}
	trikimberling(base, inclination = 0, conjugate = true){// Kimberling golden triangle
		const k = 1.3797865516812012355584834707971; // angles = {φ : φ : 1}, k = 1/(2cos(πφ/(2φ+1))
		return this.trignomon(k, base, inclination, conjugate)
	}
	trihhh(a, b, c, inclination = 0, conjugate = true){// triangle construction from altitudes (heights)
		// harmonic addition: a ÷ b = 1/(1/a + 1/b) = (a * b)/(a + b)
		function o(u, v, w){return (u*v*w)/(u*v + v*w + w*u);} // u ÷ v ÷ w = 1/(1/u + 1/v + 1/w)
		this.altitude = {a: a, b: b, c: c};
		var I = o( a,  b,  c); // inradius I
		var A = o(-a,  b,  c); // exradius A
		var B = o( a, -b,  c); // exradius B
		var C = o( a,  b, -c); // exradius C
		this.success = (A > 0) && (B > 0) && (C > 0);
		if (this.success){
			var S = 2 * sqrt(I * A * B * C); // double area
			this.trisss(S/a, S/b, S/c, inclination, conjugate); // construct from sides
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
		this.success = !(A.is(B) || C.isInf); // fail if A = B or |AB| = diameter
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
		var z = new complex(C).lens(A, B);
		var u = abs(C.r - C.zdist(z.z1));
		var v = abs(C.r - C.zdist(z.z2));
		return this.iff(z.z1, z.z2, u < v);
	}
	tripolarfun(a, b, c, f){
		return this.tripolarxiy(f(a, b, c), f(b, c, a), f(c, a, b));
	}
	toString(p = 0){// for debug
		const i = 'і';
		var z = new complex(this); if (p != 0) z.round(abs(p));
		var n = false; if (p < 0) n = z.x < 0; if (n) z.neg;
		var x = z.x.toString(), y = abs(z.y).toString(), s = '';
		y = (abs(z.y) == 1 ? '' : y + ' ') + i;
		if (z.is0) s = '0'; else
		if (z.isInf) s = '∞'; else
		if (z.isNan) s = '?'; else
		if (z.y == 0) s = x; else
		if (z.x == 0) 
			s = (z.y < 0 ? '-' : '') + y;
		else
			s = z.x + (z.y < 0 ?' - ':' + ') + y;
		if (n){
			if (z.x != 0 && z.y != 0) s = '('+s+')';
			s = '-'+s;
		}
		return s;
	}
	polyvalue(p){// evaluate
		var i = p.length, z = {}; this.obj(z).zero;
		while (0 < i--) this.zmul(z).zadd(new complex(p[i]));
		return this;
	}
	polyprime(p, q){// q = p'
		q.length = 0;
		for (var i = 1; i < p.length; i++)
			q.push(new complex(p[i]).mul(i).z);
	}
	polydiv(p, q){// p = p div q
		var m = p.hi, n = q.hi, l = p.length - q.length;
		if (l < 0) p = [0]; else {
			var r = new Array(l + 1);
			for (var k = l; k >= 0; k--){
				r[k] = {};
				var z = new complex(p[n + k]).zdiv(new complex(q[n])).obj(r[k]);
				for (var j = n; j >= 0; j--)
					new complex(q[j]).zmul(z).zbus(p[j + k]).obj(p[j + k]);
			} // p = mod, r = div
			p.length = r.length;
			for (var k = 0; k < r.length; k++)
				new complex(r[k]).obj(p[k]);
		}
	}
	polytrim(p){
		while (p.length > 0 && new complex(p[p.hi]).is0) p.length--; // leading zeroes trim
	}
	polyarg(p, ...arg){
		var n, i, a;  p.length = 0; 
		if (arg.length == 1 && Array.isArray(arg[0])) {// only one array parameter
			a = arg[0]; n = a.hi; // incremental copy array as complex polynom
			for (i = 0; i <= n; i++) p.push(new complex(a[i]).z);
		} else {// zero or more than one parameters or first parameter is not array
			a = arg; n = a.hi; // decremental collect parameters
			for (i = n; i >= 0; i--) p.push(new complex(a[i]).z);
		}
		if (this.isNum && !this.is0){// subtract this
			if (p.length == 0) p.push({}); // impossible
			new complex(p[0]).zsub(this).obj(p[0]);
		}
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
		function eps(s, t = {x: 0, y: 0}, e = 1e-52){// precision - Kahan summation
			var a = new complex(s).abs, b = new complex(t).abs, c = a - b, d = (c - a) + b;
			return (abs(c) <= e) && (abs(d) <= e);
		}
		var p = [], n, i; // polynom array, degree and index
		this.polyarg(p, ...arg); this.polytrim(p); // get and trim
		// for (i = p.hi; i >= p.lo; i--) console.log('p[',i,'] =',new complex(p[i]).toString(100000)); // (debug)
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
						w.xiy(1, 1).normaldev; // hit random point
						u.inf;  v.inf;
						while (i > 0) {
							i -= 4;
							u.asg(v);  v.asg(w);   // previous 2 iterations
							f.asg(w).polyvalue(p); // evaluate p(w)
							if (f.is0)  break;     // good luck - exact zero
							//if (eps(f)) break;     // good enough
							this.polyprime(p, q);  // q(z) = p'(z)
							g.asg(w).polyvalue(q); // evaluate p'(w)
							if (g.is0) {           // bad  luck - critical point
								w.xiy(1, 1).normaldev;
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
	get gamma(){// Γ function
		// local positive minimum at Γ(1.461632144968362341262659542325721328468)
		if (abs(this.y) > 1 || abs(this.x) > 128){// Legendre duplication formula
			//         2^(2z - 1)
			// Γ(2z) = ---------- Γ(z) Γ(z + 1/2)
			//            √ π
			var z = new complex(this).half;
			this.dec.pow2.div(sqrtpi)       // 2^(2z - 1) / sqrt(π)
				.zmul(new complex(z).gamma) // Γ(z)
				.zmul(z.add(1/2).gamma);    // Γ(z + 1/2)
		} else {// Γ(z + 1) = z Γ(z)
			var p = new complex(1); while (this.x >  1) {this.dec; p.zmul(this);} 
			var q = new complex(1); while (this.x <= 0) {q.zmul(this); this.inc;}
			// 0 < Re ≤ 1, |Im| ≤ 1
			if (this.x < 1 || this.y != 0) // no-trivial case
				if (this.isEq(1/2)) p.mul(sqrtpi); else
				q.zmul(this.polyvalue(tsgr));
			this.asg(p).zdiv(q);
		}
		return this;
	}
	get factorial(){// z! = Γ(z + 1)
		return this.inc.gamma;
	}
	get factorial2(){// z!! (improve for naturals)
		// Mathematica: Factorial2[z] // FunctionExpand
		var w = new complex(this).scl(pi).cos.dec // cos(πz) - 1
			.scl(0.1128956763223637161815488074737205358929) // ln(π/2)/4
			.exp.zmul(new complex(this.half).pow2); // (2/π)^(1/4 - cos(πz)/4) * 2^(z/2) 
		this.factorial.zmul(w); // (z/2)! * w
		return this;
	}
	beta(z){// Β(u, v) = Γ(u) Γ(v) / Γ(u + v)
		return this.asg(new complex(this).gamma.zmul(new complex(z).gamma).zdiv(this.zadd(z).gamma));
	}
	binomial(z){// this! / (z! (this - z)!)
		if (z.x != 1 || z.y != 0) if (z.x == 0 && z.y == 0 || this.is(z)) this.one; else
		this.asg(new complex(this).factorial.zdiv(new complex(z).factorial.zmul(this.zsub(z).factorial)));
		return this;
	}
	pochhammer(z){// Γ(this + z) / Γ(this)
		return this.asg(new complex(this).zadd(z).gamma.zdiv(this.gamma));
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
			if (this.sto == 0) this.del;
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

Object.defineProperty(Array.prototype, 'lo', {enumerable: false, configurable: false, get() { return 0; }});
Object.defineProperty(Array.prototype, 'hi', {enumerable: false, configurable: false, get() { return this.length - 1; }});

// Taylor Serie - Gamma Reciprocal
const tsgr = [0, 1, // N[CoefficientList[Series[1/Gamma[z], {z, 0, 30}], z], 24]
	 0.57721566490153286060651209008240243104215933593992, // γ (Euler–Mascheroni constant)
	-0.6558780715202538810770195151453904812798, // (γ² - π²/6)/2
	-0.0420026350340952355290039348754298187114,
	 0.1665386113822914895017007951021052357178,
	-0.0421977345555443367482083012891873913017,
	-9.62197152787697356211492e-03,
	 7.21894324666309954239501e-03,
	-1.16516759185906511211397e-03,
	-2.15241674114950972815730e-04,
	 1.28050282388116186153199e-04,
	-2.01348547807882386556894e-05,
	-1.25049348214267065734536e-06,
	 1.13302723198169588237413e-06,
	-2.05633841697760710345015e-07,
	 6.11609510448141581786250e-09,
	 5.00200764446922293005567e-09,
	-1.18127457048702014458813e-09,
	 1.04342671169110051049154e-10,
	 7.78226343990507125404994e-12,
	-3.69680561864220570818782e-12,
	 5.10037028745447597901548e-13,
	-2.05832605356650678322243e-14,
	-5.34812253942301798237002e-15,
	 1.22677862823826079015889e-15,
	-1.18125930169745876951376e-16,
	 1.18669225475160033257978e-18,
	 1.41238065531803178155580e-18,
	-2.29874568443537020659248e-19,
	 1.71440632192733743338396e-20,
	 1.33735173049369311486478e-22];
	 

// Не могу више да куцам Math; ја сам паскал програмер.
const min = Math.min, max = Math.max;
const sqrt = Math.sqrt, phi = (sqrt(5) + 1)/2; // φ = 1.6180339887498948482045868343656
// https://www.youtube.com/watch?v=ZPv1UV0rD8U
const pi = Math.PI, tau = 2*pi, sqrtpi = sqrt(pi); // π = 3.1415926535897932384626433832795, τ = 2π, sqrt(π)
const exp = Math.exp, log = Math.log, pow = Math.pow;
const sin = Math.sin, cos = Math.cos, tan = Math.tan;
const asin = Math.asin, acos = Math.acos;
const atan = function(s, c = 1){return s instanceof Object ? Math.atan2(s.y, s.x) : Math.atan2(s, c);}
const sqr = function(z){return z instanceof Object ? {x: z.x * z.x - z.y * z.y, y: z.x * z.y * 2} : z * z;}
const abs = function(z){return z instanceof Object ? Math.hypot(z.x, z.y) : Math.abs(z)};
const random = function(a, b){// uniform deviate in given range
	switch (arguments.length) {
		case  0: return Math.random();     // [0, 1)
		case  1: return random() * a;      // [0, a) = a * [0, 1)
		default: return random(b - a) + a; // [a, b) = a + [0, b - a)
	}
}
const xiy = function(x, y = 0){return {x: Number(x), y: Number(y)};} // shortcut
const cis = function(t, a = 1, b = a){// polar ellipse
	if (a instanceof Object) return cis(t, a.x, a.y); else
	var c = cos(t), s = sin(t), r = a == b ? a : a * b / Math.hypot(a * s, b * c);
	return {x: r * c, y: r * s};
}
const hav = function(x){return sqr(sin(x/2));};
const ahav = function(x){return 2*asin(sqrt(x));};
const sind = function(x){return sin(pi * x/180);}
const cosd = function(x){return cos(pi * x/180);}
const atand = function(s, c = 1){return atan(s, c) * 180/pi;}
const hsi2rgb = function(h, s, i){// hue in degrees, saturation and intensity = [0, 1]
	function div(x, y){return Math.floor(x/y);}
	function mod(x, y){return x - y * div(x, y);}
	function cen(x){if (x <= 0) x = 0; else if (x >= 1) x = 1; return x;} // censor x to [0, 1]
	var r, g, b, a, x, y, z;
	h = mod(h, 360); i = cen(i) * 85; s = cen(s) * i; // 255/3 = 85, s * i = chroma
	a = mod(h, 120); h = div(h, 120);
	x = i - s;
	z = i + s * cosd(a)/cosd(a - 60);
	y = 3 * i - (x + z);
	switch (h){
		case 0: b = x; g = y; r = z; break;
		case 1: b = y; g = z; r = x; break;
		case 2: b = z; g = x; r = y; break;
	}
	r = Math.round(r); g = Math.round(g); b = Math.round(b);
	return {r: r, g: g, b: b};
}


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
	this.begin.zmoveTo(z1).zlineTo(z2).stroke(); return this;
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
	this.arcTo(z1.x, z1.y, z2.x, z2.y, abs(r)); return this;
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
	this.rotate(atan(z)); return this;
}

CanvasRenderingContext2D.prototype.ztranslate = function(z){
	this.translate(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zsetTransform = function(scale, skew, move){
	this.setTransform(scale.x, skew.x, skew.y, scale.y, move.x, move.y); return this;
}


CanvasRenderingContext2D.prototype.circle = function (x, y, r) {
	this.arc(x, y, abs(r), 0, tau); return this;
}

CanvasRenderingContext2D.prototype.zcircle = function (z, r) {
	return this.circle(z.x, z.y, r);
}

CanvasRenderingContext2D.prototype.ccircle = function (c) {
	return this.zcircle(c, c.r);
}

CanvasRenderingContext2D.prototype.zellipse = 
function(z, a, b, rot = 0, s = 0, e = tau, c = false){
	this.ellipse(z.x, z.y, abs(a), abs(b), rot, s, e, c); return this;
}

CanvasRenderingContext2D.prototype.eellipse = 
function(z, s = 0, e = tau, c = false){
	return this.zellipse(z, z.a, z.b, z.o, s, e, c);
}

CanvasRenderingContext2D.prototype.fellipse = 
function(F1x, F1y, F2x, F2y, b, s = 0, e = tau, c = false){
	var dx = F1x - F2x, dy = F1y - F2y, x = F1x - dx/2, y = F1y - dy/2; o = atan(dy, dx);
	var d = dx*dx + dy*dy, a = sqrt(b*b + d); b = abs(b);
	this.ellipse(x, y, a, b, o, s, e, c); return this;
}


CanvasRenderingContext2D.prototype.superellipse = 
function(shape = 2, xpos, ypos, xradius = 1, yradius = xradius, azimuth = 0, symmetry = 4, u = shape, v = u){
	var steps = symmetry * max(9, xradius + yradius);
	if (steps > 0){
		var z = new complex(), theta;
		for (var i = 0; i < steps; i++){
			theta = i * tau / steps;
			z.superellipse(shape, theta, xradius, yradius, symmetry, u, v)
				.rotate(azimuth).add(xpos, ypos);
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
	return this.superellipse(shape, xpos, ypos, radius, radius / phi, azimuth, symmetry, u, v);
}

CanvasRenderingContext2D.prototype.zsupertable = 
function(shape = 2, z, radius = 1, azimuth = 0, symmetry = 4, u = shape, v = u){
	return this.supertable(shape, z.x, z.y, radius, azimuth, symmetry, u, v);
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

CanvasRenderingContext2D.prototype.hsia = function(h, s, i, a = 1){
	var c = hsi2rgb(h, s, i); return this.rgba(c.r, c.g, c.b, a);
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

// ... to be continued
