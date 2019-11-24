class complex {
/*
	Complex arithmetic with canvas extension
	
	JavaScript implementation by Србислав Д. Нешић, November 2019
	srbislav.nesic@gmail.com

	Euler identity e ^ (i * π) + 1 = 0: new complex(Math.PI).muli.exp.add(1);
*/
	constructor(u, v){
		this.x = 0; this.y = 0;
		if (arguments.length > 1) this.xiy(u, v); else
		if (arguments.length > 0) 
		if (typeof u === 'object') this.asg(u); else this.xiy(u);
	}
	xiy(x, y = 0){
		this.x = Number(x); this.y = Number(y); return this;
	}
	cis(rho, theta){// r * e^(ia) = r cis a = r * (cos a + i sin a)
		if (arguments.length < 2) {theta = rho; rho = 1;}
		return this.xiy(rho * Math.cos(theta), rho * Math.sin(theta));
	}
	asg(z){
		return this.xiy(z.x, z.y);
	}
	get sqrabs(){// |z|^2 = z * z'
		return this.x * this.x + this.y * this.y;
	}
	get abs(){
		if (this.x == 0 || this.y == 0)
			return Math.abs(this.x + this.y);
		else
			return Math.sqrt(this.sqrabs);
	}
	get arg(){
		return Math.atan2(this.y, this.x);
	}
	obj(z){
		z.x = this.x; z.y = this.y; return this;
	}
	isEqual(x, y = 0){
		return this.x == x && this.y == y;
	}
	isSame(z){
		return this.isEqual(z.x, z.y);
	}
	get isZero(){
		return this.isEqual(0);
	}
	get i(){
		return this.xiy(0, 1);
	}
	get real(){
		return this.xiy(this.x, 0);
	}
	get imag(){
		return this.xiy(0, this.y);
	}
	get conjg(){
		return this.xiy(this.x, -this.y);
	}
	get pos(){// 1st quadrant
		return this.xiy(Math.abs(this.x), Math.abs(this.y));
	}
	get neg(){// -z = z rotate 180
		return this.xiy(-this.x, -this.y);
	}
	get muli(){// z * i = z rotate 90
		return this.xiy(-this.y, this.x);
	}
	get divi(){// z / i = z rotate -90 = - i * z
		return this.xiy(this.y, -this.x);
	}
	rect(z){
		return (this.x - z.x) * (this.y - z.y);
	}
	trap(z){
		return (this.x - z.x) * (this.y + z.y) / 2;
	}
	adv(x, y = x){
		return this.xiy(this.x + x, this.y + y);
	}
	zadv(z){
		return this.adv(z.x, z.y);
	}
	scl(x, y = x){
		return this.xiy(this.x * x, this.y * y);
	}
	zscl(z){
		return this.scl(z.x, z.y);
	}
	add(x, y = 0){
		return this.xiy(this.x + x, this.y + y);
	}
	zadd(z){
		return this.add(z.x, z.y);
	}
	sub(x, y = 0){
		return this.xiy(this.x - x, this.y - y);
	}
	zsub(z){
		return this.sub(z.x, z.y);
	}
	mul(x, y = 0){
		return this.xiy(this.x * x - this.y * y, this.x * y + this.y * x);
	}
	zmul(z){
		return this.mul(z.x, z.y);
	}
	div(x, y = 0){
		if (y == 0) return this.xiy(this.x / x, this.y / x); else 
		if (x == 0) return this.xiy(this.x / y, this.y / y).divi; else
			return this.mul(x, -y).div(x * x + y * y);
	}
	zdiv(z){
		return this.div(z.x, z.y);
	}
	get recip(){
		if (this.y == 0) return this.xiy(1 / this.x); else 
		if (this.x == 0) return this.xiy(1 / this.y).divi; else
			return this.conjg.div(this.sqrabs);
	}
	get unit(){// unit vector
		var d = this.abs;
		return d == 0 ? this : this.div(d);
	}
	rotate(angle){
		var h = Math.PI/2;
		if (angle == 0) return this; else 
		if (angle == h) return this.muli; else 
		if (-angle == h) return this.divi; else 
		if (Math.abs(angle) == Math.PI) return this.neg; else
			return this.mul(Math.cos(angle), Math.sin(angle));
	}
	zrotate(z){
		var h = z.x * z.x + z.y * z.y;
		return h == 0 ? this : this.zmul(z).div(Math.sqrt(h));
	}
	vrotate(z1, z2){
		return this.zrotate(new complex(z2).zsub(z1));
	}
	about(x, y, angle){// rotate this about point (x, y)
		return this.sub(x, y).rotate(angle).add(x, y);
	}
	zabout(z, angle){
		return this.about(z.x, z.y, angle);
	}
	times(z1, z2 = {x: 0, y: 0}){// this = z1 + times * (z2 - z1)
		return this.zsub(z1).zdiv(new complex(z2).zsub(z1)); 
	}
	inter(z1, z2 = {x: 0, y: 0}){// inversion of times
		return this.zmul(new complex(z2).zsub(z1)).zadd(z1);
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
	isotonic(z1, z2 = {x: 0, y: 0}){
		return this.opposite(new complex(z1).halfway(z2));
	}
	perp(z1, z2 = {x: 0, y: 0}){
		return this.zadd(new complex(z2).zsub(z1).muli);
	}
	intersection(z1, z2, z3, z4){// intersection point of lines z1--z2 and z3--z4
		function cross(p, q){return p.x * q.y - p.y * q.x;}
		var u = new complex(z2).zsub(z1);
		var v = new complex(z4).zsub(z3);
		var d = cross(u, v);
		if (d == 0){// parallel
			return this.xiy(1/0, 1/0); // complex infinity
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
	zdist(z){
		return this.dist(z.x, z.y);
	}
	azimuth(x, y){// angle of this -- (0, 0) -- (x, y)
		return new complex(x, y).zdiv(this).arg;
	}
	zazimuth(z){
		return azimuth(z.x, z.y);
	}
	linedist(z1, z2 = {x: 0, y: 0}){// directed distance from line z1--z2 
		return new complex(z1).zdist(z2) * new complex(this).times(z1, z2).y;
	}
	toward(x, y, len = 1){
		if (this.isEqual(x, y))
			return this;
		else
			return this.zadd(new complex(x, y).zsub(this).unit.mul(len));
	}
	ztoward(z, len = 1){
		return this.toward(z.x, z.y, len);
	}
	oncircle(circle){
		if (this.isSame(circle))
			return this;
		else
			return this.asg(new complex(circle).ztoward(this, circle.r));
	}
	isotonicconjg(A, B, C){
		return this.intersection(
			A, new complex().intersection(A, this, B, C).isotonic(B, C), 
			B, new complex().intersection(B, this, C, A).isotonic(C, A));
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
	linecircle(point1, point2, circle){// line - circle intersection
		function cross(p, q){return p.x * q.y - p.y * q.x;};
		var u = new complex(point1).zsub(circle);
		var v = new complex(point2).zsub(circle);
		var z = new complex(v).zsub(u);
		var d = z.sqrabs, a = cross(u, v);
		var b = d * circle.r * circle.r - a * a;
		this.z1 = {x: 1/0, y: 1/0}; 
		this.z2 = {x: 1/0, y: 1/0};
		this.success = b >= 0;
		if (this.success){
			u.asg(z).mul(a).divi;
			v.asg(z).mul((z.y < 0 ? -1 : 1) * Math.sqrt(b));
			new complex(u).zadd(v).div(d).zadd(circle).obj(this.z1);
			new complex(u).zsub(v).div(d).zadd(circle).obj(this.z2);
		}
		return this;
	}
	circlecircle(circle1, circle2){// circle - circle intersection
		var z = new complex(circle2).zsub(circle1);
		var d = z.sqrabs, c = Math.sqrt(d);  z.div(c);
		var u = circle1.r, v = circle2.r; u *= u; v *= v;
		var p = d + u - v, q = 4 * d * u - p * p;
		this.z1 = {x: 1/0, y: 1/0}; 
		this.z2 = {x: 1/0, y: 1/0};
		this.success = q >= 0;
		if (this.success){
			this.xiy(p, Math.sqrt(q)).div(2 * c);
			new complex(this).zmul(z).zadd(circle1).obj(this.z1);
			new complex(this).conjg.zmul(z).zadd(circle1).obj(this.z2);
		}
		return this;
	}
	circletangent(point, circle){// tangent from point to circle
		if (true){
			var z = new complex(point).halfway(circle); z.r = z.zdist(point);
			this.circlecircle(z, circle);
		} else {
			var z = new complex(point).zsub(circle);
			var d = z.sqrabs, c = Math.sqrt(d);  z.div(c);
			var a = circle.r;
			this.z1 = {x: 1/0, y: 1/0}; 
			this.z2 = {x: 1/0, y: 1/0};
			this.success = c >= a;
			if (this.success){
				var b = Math.sqrt(d - a * a);
				this.xiy(a, b).mul(a).div(c);
				new complex(this).zmul(z).zadd(circle).obj(this.z1);
				new complex(this).conjg.zmul(z).zadd(circle).obj(this.z2);
			}
		}
		this.asg(point); this.r = this.zdist(this.z1);
		return this;
	}
	get round(){
		return this.xiy(Math.round(this.x), Math.round(this.y));
	}
	dcp(x, y = 0){// Re = dot product, Im = cross product (not commutative)
		return this.conjg.mul(x, y);
	}
	zdcp(z){
		return this.dcp(z.x, z.y);
	}
	get theo(){// spiral of Theodorus
		if (this.isZero)
			return this.xiy(1);
		else
			return this.zadd(new complex(this).unit.muli);
	}
	get sqr(){
		return this.zmul(this);
	}
	get sqrt(){
		if (this.y == 0){
			if (this.x < 0)
				return this.xiy(Math.sqrt(-this.x)).muli;
			else
				return this.xiy(Math.sqrt(this.x));
		} else
		return this.cis(Math.sqrt(this.abs), this.arg/2);
	}
	get exp(){// e^(x + iy) = e^x * e^(iy) = e^x * (cos y + i sin y) = e^x cis y
		if (this.isEqual(0, Math.PI)) return this.xiy(-1); else // to Euler
		return this.cis(Math.exp(this.x), this.y);
	}
	get log(){// ln(r cis a) = ln(r * e^(ia)) = ln r + ln(e^(ia)) = ln r + ia
		if (this.isEqual(-1)) return this.xiy(0, Math.PI); else // ditto
		return this.xiy(Math.log(this.abs), this.arg);
	}
	pow(x, y = 0){
		if (x == 0 && y == 0) return this.xiy(1); else
		if (this.isZero) return this; else
			return this.log.mul(x, y).exp;
	}
	zpow(z){
		return this.pow(z.x, z.y);
	}
	root(x, y = 0){
		return this.log.div(x, y).exp;
	}
	zroot(z){
		return this.root(z.x, z.y);
	}
	get sinh(){
		this.exp; return this.zsub(new complex(this).recip).div(2);
	}
	get cosh(){
		this.exp; return this.zadd(new complex(this).recip).div(2);
	}
	get tanh(){
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
	get sinc(){
		var z = new complex(this);
		return this.isZero ? this.xiy(1) : this.sin.zdiv(z);
	}
	get asinh(){// pricipal value
		var z = new complex(this);
		return this.sqr.add(1).sqrt.zadd(z).log;
	}
	get acosh(){// principal value
		var z = new complex(this);
		return this.sqr.sub(1).sqrt.zadd(z).log;
	}
	get atanh(){
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
	arithSpiral(a, b, t){
		return this.cis(a + b * t, t);
	}
	logSpiral(a, b, t){
		return this.cis(a * Math.exp(b * t), t);
	}
	nautilus(t){
		const a = 1, b = 0.30634896253003312211567570119977; // ln(φ) / (π/2)
		return this.logSpiral(a, b, t);
	}
	seed(n, size = 1, angle = 0){// sunflower seeds
		const f = 2.3999632297286533222315555066336; // π * (3 - sqrt(5))
		return this.cis(size * Math.sqrt(n), angle + n * f);
	}
	superellipse(shape = 2, angle, xradius = 1, yradius = xradius, symmetry = 4, u = shape, v = u){
		this.cis(angle * symmetry/4).pos.scl(1/xradius, 1/yradius);
		this.xiy(Math.pow(this.x, u), Math.pow(this.y, v));
		return this.cis(Math.pow(this.x + this.y, -1/shape), angle);
	}
	supercircle(shape = 2, angle, radius = 1, symmetry = 4, u = shape, v = u){
		return this.superellipse(shape, angle, radius, radius, symmetry, u, v);
	}
	get rectdev(){
		return this.scl(Math.random(), Math.random());
	}
	get polardev(){
		return this.zscl(
			new complex().cis(
				Math.sqrt(Math.random()), 
				2 * Math.PI * Math.random()
			)
		);
	}
	get cartdev(){// Cartesian
		return this.scl(
			2 * Math.random() - 1, 
			2 * Math.random() - 1
		);
	}
	get expdev(){
		function laplace(){
			return -Math.log(1 - Math.random()) * (Math.random() < 0.5 ? -1 : 1);
		}
		return this.scl(laplace(), laplace());
	}
	get normaldev(){
		return this.zscl(
			new complex().cis(
				Math.sqrt(-2 * Math.log(1 - Math.random())), 
				2 * Math.PI * Math.random()
			)
		);
	}
	get expgaussdev(){
		return this.zscl(
			new complex(
				-Math.log(1 - Math.random()) * (Math.random() < 0.5 ? -1 : 1),
				Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.sin(2 * Math.PI * Math.random())
			)
		);
	}
	get poissondev(){
		function poisson(lambda){
			if (lambda == 0) return 0; else {
				var l = Math.exp(-Math.abs(lambda)), r = -1, p = 1;
				do { r++; p *= Math.random(); } while (p > l);
				if (lambda < 0) r = -r;
				return r;
			};
		}
		return this.xiy(poisson(this.x), poisson(this.y));
	}
	trivertex(vertexA, vertexB, vertexC, changed = true){// triangle ABC centers
		var A = new complex(vertexA), B = new complex(vertexB), C = new complex(vertexC);
		this.vertex = {A: {x: A.x, y: A.y}, B: {x: B.x, y: B.y}, C: {x: C.x, y: C.y}};
		this.box = {
			min: {x: Math.min(A.x, B.x, C.x), y: Math.min(A.y, B.y, C.y)}, 
			max: {x: Math.max(A.x, B.x, C.x), y: Math.max(A.y, B.y, C.y)}
		}; 
		this.box.size = {x: this.box.max.x - this.box.min.x, y: this.box.max.y - this.box.min.y};
		if (changed) this.side = {a: B.zdist(C), b: C.zdist(A), c: A.zdist(B)};
		var a = this.side.a, b = this.side.b, c = this.side.c;
		this.success = (a < b + c) && (b < c + a) && (c < a + b);
		if (this.success){
			var o = new complex(A).trap(B) + new complex(B).trap(C) + new complex(C).trap(A);
			this.direction = Math.sign(o);
			var P = a + b + c,  s = P/2,  D = Math.abs(o); // Δ = D
			this.height = {a: 2*D/a, b: 2*D/b, c: 2*D/c};
			this.perimeter = P; this.area = D; this.semi = s;
			this.omega = Math.atan2(4*D, a*a + b*b + c*c); // Brocard angle
			this.I = {r: D/s}; // X1 - incircle (weighted average (Aa + Bb + Cc)/(a + b + c))
			this.xiy(0).zadd(new complex(A).mul(a)).zadd(new complex(B).mul(b)).zadd(new complex(C).mul(c)).div(P).obj(this.I);
			this.I.A = {}; this.asg(this.I).ortho(B, C).obj(this.I.A); // incircle
			this.I.B = {}; this.asg(this.I).ortho(C, A).obj(this.I.B); // contact
			this.I.C = {}; this.asg(this.I).ortho(A, B).obj(this.I.C); // points
			var p = a*b + b*c + c*a; this.I.a = this.I.r * Math.sqrt(p*p - a*b*c*s - p*s*s)/(p - s*s); // Adams radius
			this.G = {}; // X2 - centroid 
			this.xiy(0).zadd(A).zadd(B).zadd(C).div(3).obj(this.G);
			this.O = {}; // X3 - circumcircle
			this.O.r = this.bisection(A, B, C).obj(this.O).zdist(A);
			this.H = {}; // X4 - orthocenter (2OG = GH)
			this.asg(this.O).anticomplement(this.G).obj(this.H); // allways with respect to G
			this.angle = {A: Math.asin(a/2/this.O.r), B: Math.asin(b/2/this.O.r)};
			this.angle.C = Math.PI - (this.angle.A + this.angle.B);
			this.N = {}; // X5 - nine-point circle
			this.asg(this.O).halfway(this.H).obj(this.N); this.N.r = this.O.r/2;
			this.K = {}; // X6 - symmedian point (isogonal conjugate of G, Lemoine)
			this.asg(this.G).isogonalconjg(A, B, this.I).obj(this.K);
			this.Ge = {}; // X7 - Gergonne point (intersection of vertex--contact lines)
			this.intersection(A, this.I.A, B, this.I.B).obj(this.Ge);
			this.Na = {}; // X8 - Nagel point (isotonic conjugate of Ge, anticomplement of I)
			//this.asg(this.Ge).isotonicconjg(A, B, C).obj(this.Na);
			this.asg(this.I).anticomplement(this.G).obj(this.Na);
			this.M = {}; // X9 - Mittentpunkt (complement of Ge)
			this.asg(this.Ge).complement(this.G).obj(this.M);
			this.Sp = {}; // X10 - Spieker center (complement of I)
			this.asg(this.I).complement(this.G).obj(this.Sp); this.Sp.r = this.I.r/2;
			this.F = {}; // X11 - Feuerbach point (common tangent point of incircle and nine-point circle)
			this.asg(this.I).oncircle(this.N).obj(this.F);
			this.X12 = {}; // harmonic conjugate of X11 with respect to X1 and X5
			this.asg(this.F).harmonicconjg(this.I, this.N).obj(this.X12);
			this.L = {}; // X20 - de Longchamps point (Soddy line I--L)
			this.asg(this.H).opposite(this.O).obj(this.L);
			this.Euler = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Euler line (contains G, H, O, N, ...)
			this.Nagel = {z1: {x: 1/0, y: 1/0}, z2: {x: 1/0, y: 1/0}}; // Nagel line (contains G, I, Na, Sp, ...)
			if (a != b || b != c){ // no lines for equilateral triangle
				this.asg(this.G).oncircle(this.O).obj(this.Euler.z1).opposite(this.O).obj(this.Euler.z2);
				this.linecircle(this.I, this.G, this.O);
				this.asg(this.z1).obj(this.Nagel.z1); delete(this.z1);
				this.asg(this.z2).obj(this.Nagel.z2); delete(this.z2);
			}
			this.shield = {}; // shield circle
			this.shield.r = this.asg(this.G).halfway(this.H).obj(this.shield).zdist(this.G);
			var Z = Math.sqrt((a*a*a*a + b*b*b*b + c*c*c*c) - (a*a*b*b + b*b*c*c + c*c*a*a)), sq = a*a + b*b + c*c;
			this.Oe = {// Steiner circumellipse
				x: this.G.x, y: this.G.y,
				a: Math.sqrt(sq + Z*2)/3, b: Math.sqrt(sq - Z*2)/3, c: Math.sqrt(Z) * 2/3, 
				o: new complex(this.K).oncircle(this.shield).zsub(this.G).arg, F1: {}, F2: {}
			}; 
			this.Oe.e = this.Oe.c / this.Oe.a;  this.Oe.l = (sq - Z*2)/9 / this.Oe.a; // eccentricity, semi-latus rectum
			this.cis(this.Oe.c, this.Oe.o).zadd(this.Oe).obj(this.Oe.F1).opposite(this.Oe).obj(this.Oe.F2); // focii
			this.S = {}; // X99 - Steiner point (intersection of circumcircle and circumellipse)
			this.barycentricxiy(1/(b*b-c*c), 1/(c*c-a*a), 1/(a*a-b*b)).obj(this.S);
			this.asg(this.O);
		}
		return this;
	}
	triside(a, b, c, origin = {x: 0, y: 0}, inclination = 0){// triangle construction from sides
		this.side = {a: a, b: b, c: c};
		this.success = (a < b + c) && (b < c + a) && (c < a + b);
		if (this.success){
			var P = a + b + c,  s = P/2,  D = Math.sqrt(s * (s - a) * (s - b) * (s - c));  
			var A = new complex(0, 0);
			var B = new complex(c, 0);
			var C = new complex(b*b + c*c - a*a, 4*D).div(2*c).conjg;
			var O = new complex().bisection(A, B, C); // origin
			var z = new complex().cis(inclination).conjg;
			// translate to (0, 0), rotate by inclination, translate to origin
			A.zsub(O).zmul(z).zadd(origin);  
			B.zsub(O).zmul(z).zadd(origin);  
			C.zsub(O).zmul(z).zadd(origin);  
			this.trivertex(A, B, C, false);
		}
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
		var c = t.c * this.side.c/2;
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
		this.xiy((a + c * Math.cos(this.angle.B)) / Math.sin(this.angle.B), -c * this.direction);
		this.vrotate(this.vertex.B, this.vertex.A).zadd(this.I);
		return this;
	}
	barycentricxiy(a, b, c){// G = {1 : 1 : 1}
		return this.trilinearxiy(a/this.side.a, b/this.side.b, c/this.side.c);
	}
	tripolarxiy(a, b, c){// O = {R : R : R}
		var A = new complex(this.vertex.A); A.r = a;
		var B = new complex(this.vertex.B); B.r = b;
		var C = new complex(this.vertex.C); C.r = c;
		var z = new complex(C).circlecircle(A, B);
		var a = C.r - new complex(z.z1).zdist(C);
		var b = C.r - new complex(z.z2).zdist(C);
		if (Math.abs(a) < Math.abs(b)) this.asg(z.z1); else this.asg(z.z2);
		return this;
	}
}


Object.defineProperty(CanvasRenderingContext2D.prototype, 'width', {
	enumerable: false,
	configurable: false,
	get() { return this.canvas.width; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'height', {
	enumerable: false,
	configurable: false,
	get() { return this.canvas.height; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'sizeX', {
	enumerable: false,
	configurable: false,
	get() { return this.width; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'sizeY', {
	enumerable: false,
	configurable: false,
	get() { return this.height; }
});

Object.defineProperty(CanvasRenderingContext2D.prototype, 'size', {
	enumerable: false,
	configurable: false,
	get() { return {x: this.sizeX, y: this.sizeY}; }
});

CanvasRenderingContext2D.prototype.clear = function(){
	this.save();
	this.clearRect(0, 0, this.sizeX, this.sizeY);
	this.restore();
	return this;
}

CanvasRenderingContext2D.prototype.strokefill = function () {
	this.stroke(); this.fill();
	return this;
}

CanvasRenderingContext2D.prototype.startPath = function(){
	this.beginPath(); return this;
}

CanvasRenderingContext2D.prototype.endPath = function(){
	this.closePath(); return this;
}

CanvasRenderingContext2D.prototype.zmoveTo = function(z){
	this.moveTo(z.x, z.y); return this;
}

CanvasRenderingContext2D.prototype.zlineTo = function(z){
	this.lineTo(z.x, z.y); return this;
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

CanvasRenderingContext2D.prototype.zellipse = function(z, a, b, rot = 0, s = 0, e = 2*Math.PI, c = false){
	this.ellipse(z.x, z.y, a, b, rot, s, e, c); return this;
}

CanvasRenderingContext2D.prototype.eellipse = function(z, s = 0, e = 2*Math.PI, c = false){
	return this.zellipse(z, z.a, z.b, z.o, s, e, c);
}


CanvasRenderingContext2D.prototype.superellipse = 
function(shape = 2, xpos, ypos, xradius = 1, yradius = xradius, azimuth = 0, symmetry = 4, u = shape, v = u){
	var steps = symmetry * Math.max(9, xradius + yradius);
	var z = new complex();
	for (var i = 0; i < steps; i++){
		z.superellipse(shape, 2 * i * Math.PI/steps, xradius, yradius, symmetry, u, v).rotate(azimuth).add(xpos, ypos);
		if (i == 0)
			this.zmoveTo(z);
		else
			this.zlineTo(z); // can be improved with spline and less steps
	}
	if (steps > 0) this.endPath();
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

CanvasRenderingContext2D.prototype.triangle = function(t) {
	return this
	.zmoveTo(t.vertex.A)
	.zlineTo(t.vertex.B)
	.zlineTo(t.vertex.C)
	.endPath();
}

CanvasRenderingContext2D.prototype.oval = function (x, y, w, h, r = 0) {
	if (w < 2 * r) r = w / 2;
	if (h < 2 * r) r = h / 2;
	this.moveTo(x+r, y);
	this.arcTo(x+w, y,   x+w, y+h, r);
	this.arcTo(x+w, y+h, x,   y+h, r);
	this.arcTo(x,   y+h, x,   y,   r);
	this.arcTo(x,   y,   x+w, y,   r);
	return this.endPath();
}

CanvasRenderingContext2D.prototype.zoval = function (z1, z2, r = 0) {
	var z = new complex(z2).zsub(z1);
	return this.oval(z1.x, z1.y, z.x, z.y, r);
}
