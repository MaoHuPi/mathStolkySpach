/*
 * 2025 (c) MaoHuPi
 * mathStolkySpach/main.js
 * math objects and calculate methods about exploration of life
 */

/* basic */
class Test {
	constructor(id, tasks) {
		this.id = id;
		this.tasks = tasks;
	}
	run() {
		let errorCount = 0;
		for (let i = 0; i < this.tasks.length; i++) {
			if (!this.tasks[i]()) {
				console.error(`test (id ${this.id}) task (index ${i}) failed:\n${this.tasks[i].toString()}`);
				errorCount++;
			}
		}
		console.warn(`test (id ${this.id}) result: ` + (errorCount == 0 ? '%call pass' : `%c${errorCount} failed`), errorCount == 0 ? 'color: yellowgreen;' : 'color: pink;');
	};
}

/* math */
class Num {
	static quotient = Object.fromEntries([
		['Integer', 1],
		['Fraction', 2],
		['Float', 3]
	]);
	constructor() {
	}
	static getTargetType(...elements) {
		if (elements.length == 0) return this;
		return elements
			.map(element => {
				return [element.type, this.quotient[element.type.name]];
			})
			.sort((a, b) => b[1] - a[1])[0][0];
	}
	static add(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.add(getTargetType.from(a), getTargetType.from(b));
	}
	static sub(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.sub(getTargetType.from(a), getTargetType.from(b));
	}
	static mul(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.mul(getTargetType.from(a), getTargetType.from(b));
	}
	static div(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.div(getTargetType.from(a), getTargetType.from(b));
	}
	static pow(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.pow(getTargetType.from(a), getTargetType.from(b));
	}
	static mod(a, b) {
		let getTargetType = this.getTargetType(a, b);
		return getTargetType.mod(getTargetType.from(a), getTargetType.from(b));
	}
	static from(n) {
		if (typeof n == 'number') {
			if (n % 1 == 0) {
				return new Integer(n);
			} else if (Number.isFinite(n)) {
				return new Float(n);
			} else if (Number.isNaN(n)) {
				throw new Error('Num can\'t be NaN.');
			} else if (n == Infinity) {
				return new Sym(n);
			}
		} else if (n instanceof Num) {
			return n;
		}
	}
}
class Integer extends Num {
	constructor(integer) {
		super();
		this.integer = integer;
	}
	get type() {
		return Integer;
	}
	get value() {
		return this.integer;
	}
	copy() {
		return new this.constructor(this.integer);
	}
	toString() {
		return `${this.integer}`;
	}
	static add(a, b) {
		return new this(a.integer + b.integer);
	}
	static sub(a, b) {
		return new this(a.integer - b.integer);
	}
	static mul(a, b) {
		return new this(a.integer * b.integer);
	}
	static div(a, b) {
		if (a.integer % b.integer == 0) {
			return new Integer(a.integer / b.integer);
		} else {
			return new Fraction(a.integer, b.integer);
		}
	}
	static pow(a, b) {
		return new this(a.integer ** b.integer);
	}
	static mod(a, b) {
		return new this(a.integer % b.integer);
	}
	static factorization(n) {
		n = Math.abs(n.integer);
		let res = [];
		let middleI = Math.sqrt(n);
		for (let i = 1; i < middleI; i++) {
			if (n % i == 0) {
				res.push(i);
				res.push(n / i);
			}
		}
		if (middleI % 1 == 0 && n % middleI == 0) {
			res.push(middleI);
		}
		return res.sort((a, b) => a - b).map(m => new Integer(m));
	}
	static gcd(a, b) {
		a = Math.abs(a.integer);
		b = Math.abs(b.integer);
		let smallNum = Math.min(a, b);
		for (let i = smallNum; i > 0; i--) {
			if (a % i == 0 && b % i == 0) {
				return new Integer(i);
			}
		}
	}
	static from(n) {
		if (!(n instanceof Num)) {
			n = Num.from(n);
		}
		switch (n.type) {
			case Integer:
				return n.copy();
			case Fraction:
				return new this(Math.floor(n.value));
			case Float:
				return new this(Math.floor(n.value));
		}
	}
}
class Fraction extends Num {
	constructor(numerator, denominator) {
		super();
		this.numerator = numerator;
		if (denominator == 0) throw Error('Denominator of the fraction can\'t be zero.');
		this.denominator = denominator;
	}
	get type() {
		return Fraction;
	}
	get value() {
		return this.numerator / this.denominator;
	}
	tidy() {
		if (this.denominator !== 1) {
			let gcd = Integer.gcd(new Integer(this.numerator), new Integer(this.denominator));
			if (gcd.value !== 0) {
				this.numerator /= gcd.value;
				this.denominator /= gcd.value;
			}
			if (this.denominator < 0) {
				this.numerator *= -1;
				this.denominator *= -1;
			}
		}
	}
	copy() {
		return new this.constructor(this.numerator, this.denominator);
	}
	toString() {
		return `\\frac{${this.numerator}}{${this.denominator}}`;
	}
	static add(a, b) {
		return new this(
			a.numerator * b.denominator + b.numerator * a.denominator,
			a.denominator * b.denominator
		);
	}
	static sub(a, b) {
		return new this(
			a.numerator * b.denominator - b.numerator * a.denominator,
			a.denominator * b.denominator
		);
	}
	static mul(a, b) {
		return new this(
			a.numerator * b.numerator,
			a.denominator * b.denominator
		);
	}
	static div(a, b) {
		if (b.numerator == 0) throw Error('Denominator of the fraction can\'t be zero.');
		return new this(
			a.numerator * b.denominator,
			a.denominator * b.numerator
		);
	}
	static pow(a, b) {
		b.tidy();
		if (!(b.denominator == 1)) throw Error('Fraction.pow can\'t handle the power number which is not a positive integer.');
		if (b.value == 0) return new Integer(1);
		return new this(
			(b.numerator > 0 ? a.numerator : a.denominator) ** b.numerator,
			(b.numerator > 0 ? a.denominator : a.numerator) ** b.numerator
		);
	}
	static mod(a, b) {
		let times = Math.floor(a.value / b.value);
		this.sub(a, this.mul(b, new this(times, 1)));
	}
	static from(n) {
		if (!(n instanceof Num)) {
			n = Num.from(n);
		}
		switch (n.type) {
			case Integer:
				return new this(n.value, 1);
			case Fraction:
				return n.copy();
			case Float:
				if (n.value % 1 == 0) {
					return new this(n.value, 1);
				}
				let denominator = 10 ** (n.value.toString().split('.')[1].length);
				return new this(n.value * denominator, denominator);
		}
	}
}
class Float extends Num {
	constructor(float) {
		this.float = float;
	}
	get value() {
		return this.float;
	}
	copy() {
		return new this.constructor(this.float);
	}
	toString() {
		return `(${this.float})`;
	}
	static add(a, b) {
		return new Float(a.float + b.float);
	}
	static sub(a, b) {
		return new Float(a.float - b.float);
	}
	static mul(a, b) {
		return new Float(a.float * b.float);
	}
	static div(a, b) {
		return new Float(a.float / b.float);
	}
	static pow(a, b) {
		return new Float(a.float ** b.float);
	}
	static mod(a, b) {
		return new Float(a.float % b.float);
	}
	static from(n) {
		if (!(n instanceof Num)) {
			n = Num.from(n);
		}
		switch (n.type) {
			case Integer:
				return new this(n.value);
			case Fraction:
				return new this(n.value);
			case Float:
				return n.copy();
		}
	}
}
class Polynomial1V {
	constructor(coefficient) {
		this.coefficient = coefficient.map(n => n instanceof Num ? n : Num.from(n)); // coefficient[i] * (x ** i)
	}
	getValue(x) {
		this.tidy();
		return this.coefficient.map((n, i) => Num.mul(n, Num.pow(x, new Integer(i)))).reduce((s, n) => Num.add(s, n));
	}
	tidy() {
		if (this.coefficient.length > 1 && this.coefficient[this.coefficient.length - 1].value == 0) {
			let index = this.coefficient.length - 1;
			while (index > 0 && this.coefficient[index].value == 0) {
				index--;
			}
			this.coefficient.splice(index + 1, this.coefficient.length - index);
		}
		if (this.coefficient.length == 0) {
			this.coefficient = [0];
		}
	}
	copy() {
		return new this.constructor(this.coefficient);
	}
	toString() {
		this.tidy();
		return this.coefficient.map((n, i) => n.toString() + (i > 0 ? (i >= 10 ? `x^{${i}}` : `x^${i}`) : '')).toReversed().join(' + ');
	}
	get degree() {
		this.tidy();
		return this.coefficient.length - 1;
	}
	static longDivision(a, b, subStep = false) {
		a.tidy();
		if ((a.degree == 0 && a.coefficient[0].value == 0) || (a.degree < b.degree)) {
			let res = [subStep ? [] : new this.constructor([0]), a.copy()];
			res.toString = () => `(0) ... (${res[1].toString()})`;
			return res;
		} else {
			let newCoefficient = Num.div(a.coefficient[a.degree], b.coefficient[b.degree]);
			let nextDividendCoefficient = a.coefficient.map((n, i) => {
				if (i < a.degree - b.degree) {
					return n;
				} else {
					return Num.sub(n, Num.mul(b.coefficient[i - (a.degree - b.degree)], newCoefficient));
				}
			});
			let [tempQuotientCoefficient, reminder] = this.longDivision(new this(nextDividendCoefficient), b, true);
			tempQuotientCoefficient.push(newCoefficient);
			let res = ([subStep ? tempQuotientCoefficient : new this(tempQuotientCoefficient), reminder]);
			res.toString = () => `(${new this(tempQuotientCoefficient).toString()}) ... (${res[1].toString()})`;
			return res;
		}
	}
	static factorization(f) {
		// first-order factor with integer coefficients
		if (!(f.coefficient.filter(n => n.type == Integer).length == f.coefficient.length)) return;
		let aList = Integer.factorization(f.coefficient[f.degree]); // degree 1 coefficient
		let bList = Integer.factorization(f.coefficient[0]); // degree 0 coefficient
		bList = [...bList, ...bList.map(n => Integer.sub(new Integer(0), n))];
		let coefficientCompositions = aList.map(a => bList.map(b => {
			let gcd = Integer.gcd(a, b);
			return [Integer.div(a, gcd), Integer.div(b, gcd)];
		}))
			.reduce((s, n) => [...s, ...n]) // unpack
			.reduce((s, n) => { // remove repeated item
				if (!('id' in s)) {
					let m = s;
					s = { item: [m], id: [`${m[0].value},${m[1].value}`] }
				}
				let idOfN = `${n[0].value},${n[1].value}`;
				if (!s.id.includes(idOfN)) {
					s.item.push(n);
					s.id.push(idOfN);
				}
				return s;
			}).item;

		let fOfOne = f.coefficient.reduce((s, n) => Integer.add(s, n));
		if (fOfOne.value !== 0) {
			coefficientCompositions = coefficientCompositions.filter(composition => {
				return Integer.mod(fOfOne, Integer.add(composition[0], composition[1])).value == 0;
			});
		}
		let fOfNegOne = f.coefficient.map((n, i) => (i % 2 == 1 ? Integer.sub(new Integer(0), n) : n)).reduce((s, n) => Integer.add(s, n));
		if (fOfNegOne.value !== 0) {
			coefficientCompositions = coefficientCompositions.filter(composition => Integer.mod(fOfNegOne, Integer.add(composition[0], composition[1])).value == 0);
		}

		// console.log(coefficientCompositions);
		let factor1DList = coefficientCompositions
			.filter(composition => f.getValue(Integer.div(composition[1], composition[0])).value == 0)
			.map(composition => new Polynomial1V([composition[0], composition[1]]));
		factor1DList.forEach(factor => {
			f = Polynomial1V.longDivision(f, factor)[0];
		});

		let res = [...factor1DList, f];
		res.toString = () => '(' + res.map(p => p.toString()).join(') (') + ')';
		return res;
	}
}

/* main */
new Test(1, [
	// Integer
	() => (JSON.stringify(Integer.factorization(new Integer(255)).map(n => n.value)) == JSON.stringify([1, 3, 5, 15, 17, 51, 85, 255])),
	() => (Integer.gcd(new Integer(1080), new Integer(1920)).value == 120),
]).run();

let res;

res = Polynomial1V.longDivision(
	new Polynomial1V([-30, 19, 3, -5, 1]),
	new Polynomial1V([5, -4, 1])
).toString();
console.log(res);

res = Polynomial1V.longDivision(
	new Polynomial1V([-30, 19, 3, -5, 1]),
	new Polynomial1V([1])
).toString();
console.log(res);

console.log(Polynomial1V.factorization(new Polynomial1V([-30, 19, 3, -5, 1])).toString());
console.log(Polynomial1V.factorization(new Polynomial1V([14, -26, 21, -8, 1])).toString());
console.log(Polynomial1V.factorization(new Polynomial1V([12, -100, -1, 1])).toString());