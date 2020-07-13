#include "acmtoms.h"

using namespace SpecialFunctions;

void Toms::bratio(double a, double b, double x, double y,
	double* w, double* w1, int* ierr, bool log_p)
{
	bool do_swap = false;
	int n = 0, ierr1 = 0;
	double z = 0.0, a0 = 0.0, b0 = 0.0;
	double x0 = 0.0, y0 = 0.0, lambda = 0.0;
	double temp = 0.0;

	double eps = 2.0 * 0.5 * DBL_EPSILON;

	*w = log_p ? InfNaN::neginf() : 0.0;
	*w1 = log_p ? InfNaN::neginf() : 0.0;

	// safeguard, preventing infinite loops further down
	if (std::isnan(x) || std::isnan(y) || std::isnan(a) || std::isnan(b)) {
		*ierr = 9; return;
	}

	if (a < 0.0 || b < 0.0) { *ierr = 1; return; }
	if (a == 0.0 && b == 0.0) { *ierr = 2; return; }
	if (x < 0.0 || x > 1.0) { *ierr = 3; return; }
	if (y < 0.0 || y > 1.0) { *ierr = 4; return; }

	/* check that  'y == 1 - x' : */
	z = x + y - 0.5 - 0.5;
	if (fabs(z) > eps * 3.) { *ierr = 5; return; }

	bool a_lt_b = (a < b);
	*ierr = 0;
	if (x == 0.) goto L200;
	if (y == 0.) goto L210;

	if (a == 0.) goto L211;
	if (b == 0.) goto L201;

	eps = std::max(eps, 1e-15);
	
	if ((a_lt_b ? b : a) < eps * .001) {
		if (log_p) {
			if (a_lt_b) {
				*w = log1p(-a / (a + b)); // notably if a << b
				*w1 = log(a / (a + b));
			}
			else { // b <= a
				*w = log(b / (a + b));
				*w1 = log1p(-b / (a + b));
			}
		}
		else {
			*w = b / (a + b);
			*w1 = a / (a + b);
		}
		return;
	}

#define SET_0_noswap \
    a0 = a;  x0 = x; \
    b0 = b;  y0 = y;

#define SET_0_swap   \
    a0 = b;  x0 = y; \
    b0 = a;  y0 = x;

	if (std::min(a, b) <= 1.) { /*------------------------ a <= 1  or  b <= 1 ---- */

		do_swap = (x > 0.5);
		if (do_swap) {
			SET_0_swap;
		}
		else {
			SET_0_noswap;
		}
		/* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */
		if (b0 < std::min(eps, eps * a0)) { /* L80: */
			*w = fpser(a0, b0, x0, eps, log_p);
			temp = (*w) > -M_LN2 ? log(-rexpm1(*w)) : log1p(-exp(*w));
			*w1 = log_p ? temp : 0.5 - *w + 0.5;
			goto L_end;
		}

		if (a0 < std::min(eps, eps * b0) && b0 * x0 <= 1.) { /* L90: */
			*w1 = apser(a0, b0, x0, eps);
			goto L_end_from_w1;
		}

		bool did_bup = false;
		if (std::max(a0, b0) > 1.) { /* L20:  min(a,b) <= 1 < max(a,b)  */
			if (b0 <= 1.) goto L_w_bpser;

			if (x0 >= 0.29) /* was 0.3, PR#13786 */	goto L_w1_bpser;

			if (x0 < 0.1 && pow(x0 * b0, a0) <= 0.7)	goto L_w_bpser;

			if (b0 > 15.) {
				*w1 = 0.;
				goto L131;
			}
		}
		else {
			if (a0 >= std::min(0.2, b0))	goto L_w_bpser;

			if (pow(x0, a0) <= 0.9) 	goto L_w_bpser;

			if (x0 >= 0.3)		goto L_w1_bpser;
		}
		n = 20; /* goto L130; */
		*w1 = bup(b0, a0, y0, x0, n, eps, false); did_bup = true;
		b0 += n;
	L131:
		bgrat(b0, a0, y0, x0, w1, 15 * eps, &ierr1, false);

		if (*w1 == 0 || (0 < *w1 && *w1 < DBL_MIN)) {
			if (did_bup) { // re-do that part on log scale:
				*w1 = bup(b0 - n, a0, y0, x0, n, eps, true);
			}
			else *w1 = InfNaN::neginf(); // = 0 on log-scale
			bgrat(b0, a0, y0, x0, w1, 15 * eps, &ierr1, true);
			if (ierr1) *ierr = 10 + ierr1;
			goto L_end_from_w1_log;
		}
		
		if (ierr1) *ierr = 10 + ierr1;
		if (*w1 < 0)
			std::cout << "Warning: bratio(a = " << a << ", b = " << b
				<< ", x = " << x << "): bgrat() -> w1 = " << *w1 << std::endl;
		goto L_end_from_w1;
	}
	else {
		lambda = std::isfinite(a + b)
			? ((a > b) ? (a + b) * y - b : a - (a + b) * x)
			: a * y - b * x;
		do_swap = (lambda < 0.);
		if (do_swap) {
			lambda = -lambda;
			SET_0_swap;
		}
		else {
			SET_0_noswap;
		}

		if (b0 < 40.) {
			if (b0 * x0 <= 0.7
				|| (log_p && lambda > 650.)) // << added 2010-03; svn r51327
				goto L_w_bpser;
			else
				goto L140;
		}
		else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
			if (b0 <= 100. || lambda > b0 * 0.03)
				goto L_bfrac;

		}
		else if (a0 <= 100.) {
			goto L_bfrac;
		}
		else if (lambda > a0 * 0.03) {
			goto L_bfrac;
		}

		/* else if none of the above    L180: */
		*w = basym(a0, b0, lambda, eps * 100.0, log_p);
		temp = (*w) > -M_LN2 ? log(-rexpm1(*w)) : log1p(-exp(*w));
		*w1 = log_p ? temp : 0.5 - *w + 0.5;
		goto L_end;

	} /* else: a, b > 1 */

/* EVALUATION OF THE APPROPRIATE ALGORITHM */

L_w_bpser: // was L100
	*w = bpser(a0, b0, x0, eps, log_p);
	temp = (*w) > -M_LN2 ? log(-rexpm1(*w)) : log1p(-exp(*w));
	*w1 = log_p ? temp : 0.5 - *w + 0.5;
	goto L_end;

L_w1_bpser:  // was L110
	*w1 = bpser(b0, a0, y0, eps, log_p);
	temp = (*w1) > -M_LN2 ? log(-rexpm1(*w1)) : log1p(-exp(*w1));
	*w = log_p ? temp : 0.5 - *w1 + 0.5;
	goto L_end;

L_bfrac:
	*w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
	temp = (*w) > -M_LN2 ? log(-rexpm1(*w)) : log1p(-exp(*w));
	*w1 = log_p ? temp : 0.5 - *w + 0.5;
	goto L_end;

L140:
	/* b0 := fractional_part( b0 )  in (0, 1]  */
	n = (int)b0;
	b0 -= n;
	if (b0 == 0.) {
		--n; b0 = 1.;
	}

	*w = bup(b0, a0, y0, x0, n, eps, false);

	if (*w < DBL_MIN && log_p) {
		goto L_w_bpser;
	}
	if (x0 <= 0.7) {
		*w += bpser(a0, b0, x0, eps, /* log_p = */ false);
		goto L_end_from_w;
	}
	/* L150: */
	if (a0 <= 15.) {
		n = 20;
		*w += bup(a0, b0, x0, y0, n, eps, false);
		a0 += n;
	}

	bgrat(a0, b0, x0, y0, w, 15 * eps, &ierr1, false);
	if (ierr1) *ierr = 10 + ierr1;
	goto L_end_from_w;


	/* TERMINATION OF THE PROCEDURE */

L200:
	if (a == 0.) { *ierr = 6;    return; }
	// else:
L201:
	*w = log_p ? InfNaN::neginf() : 0.0;
	*w1 = log_p ? 0.0 : 1.0;
	return;

L210:
	if (b == 0.) { *ierr = 7;    return; }
	// else:
L211:
	*w = log_p ? 0.0 : 1.0;
	*w1 = log_p ? InfNaN::neginf() : 0.0;
	return;

L_end_from_w:
	if (log_p) {
		*w1 = log1p(-*w);
		*w = log(*w);
	}
	else {
		*w1 = 0.5 - *w + 0.5;
	}
	goto L_end;

L_end_from_w1:
	if (log_p) {
		*w = log1p(-*w1);
		*w1 = log(*w1);
	}
	else {
		*w = 0.5 - *w1 + 0.5;
	}
	goto L_end;

L_end_from_w1_log:
	if (log_p) {
		*w = (*w1) > -M_LN2 ? log(-rexpm1(*w1)) : log1p(-exp(*w1));
	}
	else {
		*w = -expm1(*w1);
		*w1 = exp(*w1);
	}
	goto L_end;


L_end:
	if (do_swap) { /* swap */
		double t = *w; *w = *w1; *w1 = t;
	}
	return;

} // bratio

double Toms::gsumln(double a, double b)
{
	double x = a + b - 2.;/* in [0, 2] */

	if (x <= 0.25)
		return gamln1(x + 1.);

	if (x <= 1.25)
		return gamln1(x) + alnrel(x);

	return gamln1(x - 1.) + log(x * (x + 1.));
} /* gsumln */

double Toms::bcorr(double a0, double b0)
{
	constexpr double c0 = 0.0833333333333333;
	constexpr double c1 = -0.00277777777760991;
	constexpr double c2 = 7.9365066682539e-4;
	constexpr double c3 = -5.9520293135187e-4;
	constexpr double c4 = 8.37308034031215e-4;
	constexpr double c5 = -0.00165322962780713;

	/* System generated locals */
	double ret_val = 0.0, r1 = 0.0;

	double a = std::min(a0, b0);
	double b = std::max(a0, b0);

	double h = a / b;
	double c = h / (h + 1.);
	double x = 1. / (h + 1.);
	double x2 = x * x;

	double s3 = x + x2 + 1.;
	double s5 = x + x2 * s3 + 1.;
	double s7 = x + x2 * s5 + 1.;
	double s9 = x + x2 * s7 + 1.;
	double s11 = x + x2 * s9 + 1.;

	r1 = 1. / b;
	double t = r1 * r1;
	double w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) *
		t + c2 * s5) * t + c1 * s3) * t + c0;
	w *= c / b;

	r1 = 1.0 / a;
	t = r1 * r1;
	ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a + w;
	return ret_val;
} /* bcorr */

double Toms::algdiv(double a, double b)
{
	constexpr double c0 = 0.0833333333333333;
	constexpr double c1 = -0.00277777777760991;
	constexpr double c2 = 7.9365066682539e-4;
	constexpr double c3 = -5.9520293135187e-4;
	constexpr double c4 = 8.37308034031215e-4;
	constexpr double c5 = -0.00165322962780713;

	double c = 0.0, d = 0.0, h = 0.0, x = 0.0;

	if (a > b) {
		h = b / a;
		c = 1. / (h + 1.0);
		x = h / (h + 1.0);
		d = a + (b - 0.5);
	}
	else {
		h = a / b;
		c = h / (h + 1.0);
		x = 1.0 / (h + 1.0);
		d = b + (a - 0.5);
	}

	double x2 = x * x;
	double s3 = x + x2 + 1.;
	double s5 = x + x2 * s3 + 1.;
	double s7 = x + x2 * s5 + 1.;
	double s9 = x + x2 * s7 + 1.;
	double s11 = x + x2 * s9 + 1.;

	double t = 1. / (b * b);
	double w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) *
		t + c2 * s5) * t + c1 * s3) * t + c0;
	w *= c / b;

	double u = d * alnrel(a / b);
	double v = a * (log(b) - 1.);
	if (u > v)
		return w - v - u;
	else
		return w - u - v;
} /* algdiv */

double Toms::gamln(double a)
{
	constexpr double d = 0.418938533204673;/* d == 0.5*(LN(2*PI) - 1) */

	constexpr double c0 = 0.0833333333333333;
	constexpr double c1 = -0.00277777777760991;
	constexpr double c2 = 7.9365066682539e-4;
	constexpr double c3 = -5.9520293135187e-4;
	constexpr double c4 = 8.37308034031215e-4;
	constexpr double c5 = -0.00165322962780713;

	if (a <= 0.8)
		return gamln1(a) - log(a); /* ln(G(a+1)) - ln(a) == ln(G(a+1)/a) = ln(G(a)) */
	else if (a <= 2.25)
		return gamln1(a - 0.5 - 0.5);

	else if (a < 10.) {
		int n = (int)(a - 1.25);
		double t = a;
		double w = 1.;
		for (int i = 1; i <= n; ++i) {
			t += -1.;
			w *= t;
		}
		return gamln1(t - 1.) + log(w);
	}
	else { /* a >= 10 */
		double t = 1. / (a * a);
		double w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
		return d + w + (a - 0.5) * (log(a) - 1.);
	}
} /* gamln */


double Toms::betaln(double a0, double b0)
{
	int n = 0;
	constexpr double e = 0.918938533204673;/* e == 0.5*LN(2*PI) */
	double a = std::min(a0, b0);
	double b = std::max(a0, b0);

	if (a < 8.) {
		if (a < 1.) {
			if (b < 8.)
				return gamln(a) + (gamln(b) - gamln(a + b));
			else
				return gamln(a) + algdiv(a, b);
		}
		double w = 0.0;
		if (a < 2.) {
			if (b <= 2.) {
				return gamln(a) + gamln(b) - gsumln(a, b);
			}

			if (b < 8.) {
				w = 0.;
				goto L40;
			}
			return gamln(a) + algdiv(a, b);
		}

		if (b <= 1e3) {
			n = (int)(a - 1.);
			w = 1.0;
			for (int i = 1; i <= n; ++i) {
				a += -1.;
				double h = a / b;
				w *= h / (h + 1.);
			}
			w = log(w);

			if (b >= 8.)
				return w + gamln(a) + algdiv(a, b);
			// else
		L40:
			// 	1 < A <= B < 8 :  reduction of B
			n = (int)(b - 1.);
			double z = 1.;
			for (int i = 1; i <= n; ++i) {
				b += -1.;
				z *= b / (a + b);
			}
			return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
		}
		else { // L50:	reduction of A when  B > 1000
			int n = (int)(a - 1.);
			w = 1.;
			for (int i = 1; i <= n; ++i) {
				a += -1.;
				w *= a / (a / b + 1.);
			}
			return log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
		}

	}
	else {
		double w = bcorr(a, b);
		double h = a / b;
		double u = -(a - 0.5) * log(h / (h + 1.0));
		double v = b * alnrel(h);
		if (u > v)
			return log(b) * -0.5 + e + w - v - u;
		else
			return log(b) * -0.5 + e + w - u - v;
	}

} /* betaln */

double Toms::psi(double x)
{
	const double piov4 = .785398163397448;
	const double dx0 = 1.461632144968362341262659542325721325;

	const double p1[7] = { .0089538502298197,4.77762828042627,
		142.441585084029,1186.45200713425,3633.51846806499,
		4138.10161269013,1305.60269827897 };
	const double q1[6] = { 44.8452573429826,520.752771467162,
		2210.0079924783,3641.27349079381,1908.310765963,
		6.91091682714533e-6 };

	const double p2[4] = { -2.12940445131011,-7.01677227766759,
		-4.48616543918019,-.648157123766197 };
	const double q2[4] = { 32.2703493791143,89.2920700481861,
		54.6117738103215,7.77788548522962 };

	int m = 0, n = 0, nq = 0;
	double w = 0.0, z = 0.0;
	double den = 0.0, sgn = 0.0, xmx0 = 0.0, upper = 0.0;

	double xmax1 = (double)INT_MAX;
	double d2 = 0.5 / 0.5 * DBL_EPSILON; /*= 0.5 / (0.5 * DBL_EPS) = 1/DBL_EPSILON = 2^52 */
	if (xmax1 > d2) xmax1 = d2;

	double xsmall = 1e-9;
	double aug = 0.0;
	if (x < 0.5) {
		if (fabs(x) <= xsmall) {

			if (x == 0.) {
				goto L_err;
			}
			aug = -1. / x;
		}
		else {
			w = -x;
			sgn = piov4;
			if (w <= 0.) {
				w = -w;
				sgn = -sgn;
			}

			if (w >= xmax1) {
				goto L_err;
			}
			nq = (int)w;
			w -= (double)nq;
			nq = (int)(w * 4.);
			w = (w - (double)nq * 0.25) * 4.;

			n = nq / 2;
			if (n + n != nq) {
				w = 1. - w;
			}
			z = piov4 * w;
			m = n / 2;
			if (m + m != n) {
				sgn = -sgn;
			}

			n = (nq + 1) / 2;
			m = n / 2;
			m += m;
			if (m == n) {
				if (z == 0.) {
					goto L_err;
				}
				aug = sgn * (cos(z) / sin(z) * 4.);
			}
			else { /* L140: */
				aug = sgn * (sin(z) / cos(z) * 4.);
			}
		}
		x = 1.0 - x;
	}

	if (x <= 3.) {
		den = x;
		upper = p1[0] * x;

		for (int i = 1; i <= 5; ++i) {
			den = (den + q1[i - 1]) * x;
			upper = (upper + p1[i]) * x;
		}

		den = (upper + p1[6]) / (den + q1[5]);
		xmx0 = x - dx0;
		return den * xmx0 + aug;
	}
	if (x < xmax1) {
		w = 1.0 / (x * x);
		den = w;
		upper = p2[0] * w;

		for (int i = 1; i <= 3; ++i) {
			den = (den + q2[i - 1]) * w;
			upper = (upper + p2[i]) * w;
		}

		aug = upper / (den + q2[3]) - 0.5 / x + aug;
	}
	return aug + log(x);

L_err:
	return 0.;
} /* psi */

double Toms::gamln1(double a)
{
	double w = 0.0;
	if (a < 0.6) {
		constexpr double p0 = 0.577215664901533;
		constexpr double p1 = 0.844203922187225;
		constexpr double p2 = -0.168860593646662;
		constexpr double p3 = -0.780427615533591;
		constexpr double p4 = -0.402055799310489;
		constexpr double p5 = -0.0673562214325671;
		constexpr double p6 = -0.00271935708322958;
		constexpr double q1 = 2.88743195473681;
		constexpr double q2 = 3.12755088914843;
		constexpr double q3 = 1.56875193295039;
		constexpr double q4 = 0.361951990101499;
		constexpr double q5 = 0.0325038868253937;
		constexpr double q6 = 6.67465618796164e-4;
		w = ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a + p0) /
			((((((q6 * a + q5) * a + q4) * a + q3) * a + q2) * a + q1) * a + 1.0);
		return -(a)*w;
	}
	else { /* 0.6 <= a <= 1.25 */
		constexpr double r0 = 0.422784335098467;
		constexpr double r1 = 0.848044614534529;
		constexpr double r2 = 0.565221050691933;
		constexpr double r3 = 0.156513060486551;
		constexpr double r4 = 0.017050248402265;
		constexpr double r5 = 4.97958207639485e-4;
		constexpr double s1 = 1.24313399877507;
		constexpr double s2 = 0.548042109832463;
		constexpr double s3 = 0.10155218743983;
		constexpr double s4 = 0.00713309612391;
		constexpr double s5 = 1.16165475989616e-4;
		double x = a - 0.5 - 0.5;
		w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
			(((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0);
		return x * w;
	}
} /* gamln1 */

double Toms::gam1(double a)
{
	double w = 0.0, bot = 0.0, top = 0.0;

	double t = a;
	double d = a - 0.5;

	if (d > 0.)
		t = d - 0.5;
	if (t < 0.) { /* L30: */
		const double r[9] = { -.422784335098468,-.771330383816272,
					 -.244757765222226,.118378989872749,9.30357293360349e-4,
					 -.0118290993445146,.00223047661158249,2.66505979058923e-4,
					 -1.32674909766242e-4 };
		const double s1 = .273076135303957;
		const double s2 = .0559398236957378;

		top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]
			) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
		bot = (s2 * t + s1) * t + 1.;
		w = top / bot;

		if (d > 0.)
			return t * w / a;
		else
			return a * (w + 0.5 + 0.5);

	}
	else if (t == 0) { // L10: a in {0, 1}
		return 0.;
	}
	else { /* t > 0;  L20: */
		const double p[7] = { .577215664901533,-.409078193005776,
				 -.230975380857675,.0597275330452234,.0076696818164949,
				 -.00514889771323592,5.89597428611429e-4 };
		const double q[5] = { 1.,.427569613095214,.158451672430138,
				 .0261132021441447,.00423244297896961 };

		top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]
			) * t + p[1]) * t + p[0];
		bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
		w = top / bot;

		if (d > 0.) /* L21: */
			return t / a * (w - 0.5 - 0.5);
		else
			return a * w;
	}
} /* gam1 */

double Toms::erfc1(int ind, double x)
{
	const double c = 0.564189583547756;
	const double a[5] = { 7.7105849500132e-5,-.00133733772997339,
		.0323076579225834,.0479137145607681,.128379167095513 };
	const double b[3] = { .00301048631703895,.0538971687740286,
		.375795757275549 };
	const double p[8] = { -1.36864857382717e-7,.564195517478974,
		7.21175825088309,43.1622272220567,152.98928504694,
		339.320816734344,451.918953711873,300.459261020162 };
	const double q[8] = { 1.,12.7827273196294,77.0001529352295,
		277.585444743988,638.980264465631,931.35409485061,
		790.950925327898,300.459260956983 };
	const double r[5] = { 2.10144126479064,26.2370141675169,
		21.3688200555087,4.6580782871847,.282094791773523 };
	const double s[4] = { 94.153775055546,187.11481179959,
		99.0191814623914,18.0124575948747 };

	double ret_val = 0.0;
	double e = 0.0, t = 0.0, w = 0.0, bot = 0.0, top = 0.0;

	double ax = fabs(x);

	if (ax <= 0.5) {
		double t = x * x,
			top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.,
			bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
		ret_val = 0.5 - x * (top / bot) + 0.5;
		if (ind != 0) {
			ret_val = exp(t) * ret_val;
		}
		return ret_val;
	}

	if (ax <= 4.) {
		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
			+ p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
			+ q[5]) * ax + q[6]) * ax + q[7];
		ret_val = top / bot;

	}
	else {
		if (x <= -5.6) {
			ret_val = 2.;
			if (ind != 0) {
				ret_val = exp(x * x) * 2.;
			}
			return ret_val;
		}
		if (ind == 0 && (x > 100. || x * x > -exparg(1))) {
			return 0.;
		}

		t = 1. / (x * x);
		top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
		bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
		ret_val = (c - t * top / bot) / ax;
	}

	if (ind != 0) {
		if (x < 0.)
			ret_val = exp(x * x) * 2. - ret_val;
	}
	else {
		w = x * x;
		t = w;
		e = w - t;
		ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
		if (x < 0.)
			ret_val = 2. - ret_val;
	}
	return ret_val;

} /* erfc1 */

double Toms::erf__(double x)
{
	const double c = 0.564189583547756;
	const double a[5] = { 7.7105849500132e-5,-.00133733772997339,
		.0323076579225834,.0479137145607681,.128379167095513 };
	const double b[3] = { .00301048631703895,.0538971687740286,
		.375795757275549 };
	const double p[8] = { -1.36864857382717e-7,.564195517478974,
		7.21175825088309,43.1622272220567,152.98928504694,
		339.320816734344,451.918953711873,300.459261020162 };
	const double q[8] = { 1.,12.7827273196294,77.0001529352295,
		277.585444743988,638.980264465631,931.35409485061,
		790.950925327898,300.459260956983 };
	const double r[5] = { 2.10144126479064,26.2370141675169,
		21.3688200555087,4.6580782871847,.282094791773523 };
	const double s[4] = { 94.153775055546,187.11481179959,
		99.0191814623914,18.0124575948747 };

	/* Local variables */
	double t = 0.0, x2 = 0.0, bot = 0.0, top = 0.0;

	double ax = fabs(x);
	if (ax <= 0.5) {
		t = x * x;
		top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
		bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;

		return x * (top / bot);
	}

	if (ax <= 4.) { //  |x| in (0.5, 4]
		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
			+ p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
			+ q[5]) * ax + q[6]) * ax + q[7];
		double R = 0.5 - exp(-x * x) * top / bot + 0.5;
		return (x < 0) ? -R : R;
	}

	if (ax >= 5.8) {
		return x > 0 ? 1 : -1;
	}

	x2 = x * x;
	t = 1. / x2;
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
	t = (c - top / (x2 * bot)) / ax;
	double R = 0.5 - exp(-x2) * t + 0.5;
	return (x < 0) ? -R : R;
} /* erf */

double Toms::rlog1(double x)
{
	const double a = 0.0566749439387324;
	const double b = 0.0456512608815524;
	const double p0 = 0.333333333333333;
	const double p1 = -0.224696413112536;
	const double p2 = 0.00620886815375787;
	const double q1 = -1.27408923933623;
	const double q2 = 0.354508718369557;

	double h = 0.0, r = 0.0, t = 0.0, w = 0.0, w1 = 0.0;
	if (x < -0.39 || x > 0.57) { /* direct evaluation */
		w = x + 0.5 + 0.5;
		return x - log(w);
	}
	/* else */
	if (x < -0.18) { /* L10: */
		h = x + .3;
		h /= .7;
		w1 = a - h * .3;
	}
	else if (x > 0.18) { /* L20: */
		h = x * .75 - .25;
		w1 = b + h / 3.;
	}
	else { /*		Argument Reduction */
		h = x;
		w1 = 0.;
	}

	r = h / (h + 2.);
	t = r * r;
	w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
	return t * 2. * (1. / (1. - r) - r * w) + w1;

} /* rlog1 */

double Toms::alnrel(double a)
{
	if (fabs(a) > 0.375)
		return log(1. + a);

	const double p1 = -1.29418923021993;
	const double p2 = .405303492862024;
	const double p3 = -.0178874546012214;
	const double q1 = -1.62752256355323;
	const double q2 = .747811014037616;
	const double q3 = -.0845104217945565;
	
	double t = a / (a + 2.);
	double t2 = t * t;
	double w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) /
		(((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
	return t * 2. * w;
} /* alnrel */

double Toms::rexpm1(double x)
{
	const double p1 = 9.14041914819518e-10;
	const double p2 = .0238082361044469;
	const double q1 = -.499999999085958;
	const double q2 = .107141568980644;
	const double q3 = -.0119041179760821;
	const double q4 = 5.95130811860248e-4;

	if (fabs(x) <= 0.15) {
		return x * (((p2 * x + p1) * x + 1.) /
			((((q4 * x + q3) * x + q2) * x + q1) * x + 1.));
	}
	else { /* |x| > 0.15 : */
		double w = exp(x);
		if (x > 0.)
			return w * (0.5 - 1. / w + 0.5);
		else
			return w - 0.5 - 0.5;
	}
} /* rexpm1 */

double Toms::esum(int mu, double x, bool give_log)
{
	if (give_log)
		return x + (double)mu;

	double w = 0.0;
	if (x > 0.) { /* L10: */
		if (mu > 0)  return exp((double)mu) * exp(x);
		w = mu + x;
		if (w < 0.) return exp((double)mu) * exp(x);
	}
	else { /* x <= 0 */
		if (mu < 0)  return exp((double)mu) * exp(x);
		w = mu + x;
		if (w > 0.) return exp((double)mu) * exp(x);
	}
	return exp(w);
} /* esum */

double Toms::exparg(int l)
{
	const double lnb = 0.69314718055995;
	int m = (l == 0) ? DBL_MAX_EXP : DBL_MIN_EXP - 1;

	return m * lnb * 0.99999;
} /* exparg */

double Toms::basym(double a, double b, double lambda, double eps, bool log_p)
{
	constexpr auto num_IT = 20;

	const double e0 = 1.12837916709551;/* e0 == 2/sqrt(pi) */
	const double e1 = 0.353553390593274;/* e1 == 2^(-3/2)   */
	const double ln_e0 = 0.120782237635245; /* == ln(e0) */

	double a0[num_IT + 1], b0[num_IT + 1], c[num_IT + 1], d[num_IT + 1];

	double f = a * rlog1(-lambda / a) + b * rlog1(lambda / b);
	double t = 0.0;
	if (log_p)
		t = -f;
	else {
		t = exp(-f);
		if (t == 0.) {
			return 0; /* once underflow, always underflow .. */
		}
	}
	double z0 = sqrt(f),
		z = z0 / e1 * 0.5,
		z2 = f + f,
		h, r0, r1, w0;

	if (a < b) {
		h = a / b;
		r0 = 1. / (h + 1.);
		r1 = (b - a) / b;
		w0 = 1. / sqrt(a * (h + 1.));
	}
	else {
		h = b / a;
		r0 = 1. / (h + 1.);
		r1 = (b - a) / a;
		w0 = 1. / sqrt(b * (h + 1.));
	}

	a0[0] = r1 * 0.66666666666666663;
	c[0] = a0[0] * -0.5;
	d[0] = -c[0];
	double j0 = 0.5 / e0 * erfc1(1, z0),
		j1 = e1,
		sum = j0 + d[0] * w0 * j1;

	double s = 1.,
		h2 = h * h,
		hn = 1.,
		w = w0,
		znm1 = z,
		zn = z2;
	for (int n = 2; n <= num_IT; n += 2) {
		hn *= h2;
		a0[n - 1] = r0 * 2. * (h * hn + 1.) / (n + 2.);
		int np1 = n + 1;
		s += hn;
		a0[np1 - 1] = r1 * 2. * s / (n + 3.);

		for (int i = n; i <= np1; ++i) {
			double r = (i + 1.) * -0.5;
			b0[0] = r * a0[0];
			for (int m = 2; m <= i; ++m) {
				double bsum = 0.;
				for (int j = 1; j <= m - 1; ++j) {
					int mmj = m - j;
					bsum += (j * r - mmj) * a0[j - 1] * b0[mmj - 1];
				}
				b0[m - 1] = r * a0[m - 1] + bsum / m;
			}
			c[i - 1] = b0[i - 1] / (i + 1.);

			double dsum = 0.;
			for (int j = 1; j <= i - 1; ++j) {
				dsum += d[i - j - 1] * c[j - 1];
			}
			d[i - 1] = -(dsum + c[i - 1]);
		}

		j0 = e1 * znm1 + (n - 1.) * j0;
		j1 = e1 * zn + n * j1;
		znm1 = z2 * znm1;
		zn = z2 * zn;
		w *= w0;
		double t0 = d[n - 1] * w * j0;
		w *= w0;
		double t1 = d[np1 - 1] * w * j1;
		sum += t0 + t1;
		if (fabs(t0) + fabs(t1) <= eps * sum) {
			break;
		}
	}

	if (log_p)
		return ln_e0 + t - bcorr(a, b) + log(sum);
	else {
		double u = exp(-bcorr(a, b));
		return e0 * t * u * sum;
	}

} /* basym_ */

double Toms::grat_r(double a, double x, double log_r, double eps)
{
	if (a * x == 0.) { /* L130: */
		if (x <= a) {
			/* L100: */ return exp(-log_r);
		}
		else {
			/* L110:*/  return 0.;
		}
	}
	else if (a == 0.5) { // e.g. when called from pt()
	/* L120: */
		if (x < 0.25) {
			double p = erf__(sqrt(x));
			return (0.5 - p + 0.5) * exp(-log_r);

		}
		else { // 2013-02-27: improvement for "large" x: direct computation of q/r:
			double sx = sqrt(x);
			double q_r = erfc1(1, sx) / sx * M_SQRT_PI;
			return q_r;
		}
	}
	else if (x < 1.1) { /* L10:  Taylor series for  P(a,x)/x^a */
		double an = 3.0;
		double c = x;
		double sum = x / (a + 3.0);
		double tol = eps * 0.1 / (a + 1.0), t;
		do {
			an += 1.;
			c *= -(x / an);
			t = c / (a + an);
			sum += t;
		} while (fabs(t) > tol);

		double j = a * x * ((sum / 6. - 0.5 / (a + 2.)) * x + 1. / (a + 1.)),
			z = a * log(x),
			h = gam1(a),
			g = h + 1.;

		if ((x >= 0.25 && (a < x / 2.59)) || (z > -0.13394)) {
			// L40:
			double l = rexpm1(z),
				q = ((l + 0.5 + 0.5) * j - l) * g - h;
			if (q <= 0.) {
				/* L110:*/ return 0.;
			}
			else {
				return q * exp(-log_r);
			}
		}
		else {
			double p = exp(z) * g * (0.5 - j + 0.5);
			return /* q/r = */ (0.5 - p + 0.5) * exp(-log_r);
		}
	}
	else {
		double a2n_1 = 1.,
			a2n = 1.,
			b2n_1 = x,
			b2n = x + (1. - a),
			c = 1., am0, an0;

		do {
			a2n_1 = x * a2n + c * a2n_1;
			b2n_1 = x * b2n + c * b2n_1;
			am0 = a2n_1 / b2n_1;
			c += 1.;
			double c_a = c - a;
			a2n = a2n_1 + c_a * a2n;
			b2n = b2n_1 + c_a * b2n;
			an0 = a2n / b2n;
		} while (fabs(an0 - am0) >= eps * an0);

		return /* q/r = (r * an0)/r = */ an0;
	}
} /* grat_r */

void Toms::bgrat(double a, double b, double x, double y, double* w,
	double eps, int* ierr, bool log_w)
{
	constexpr auto n_terms_bgrat = 30;
	double c[n_terms_bgrat], d[n_terms_bgrat];
	double bm1 = b - 0.5 - 0.5;
	double nu = a + bm1 * 0.5;
	double lnx = (y > 0.375) ? log(x) : alnrel(-y);
	double z = -nu * lnx; // z =: u in (9.1) of D.&M.(1992)

	if (b * z == 0.) {
		std::cout << "Warning: bgrat(a = " << a << ", b = " << b << ", x = "
			<< x << ", y = " << y << "): z = " << z << " == 0 underflow, hence "
			<< "inaccurate pbeta()" << std::endl;
		*ierr = 1;
		return;
	}

	double log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx;
	double log_u = log_r - (algdiv(b, a) + b * log(nu));
	double u = exp(log_u);

	if (log_u == InfNaN::neginf()) {
		*ierr = 2;
		return;
	}

	bool u_0 = (u == 0.); // underflow --> do work with log(u) == log_u !
	double l = // := *w/u .. but with care: such that it also works when u underflows to 0:
		log_w
		? ((*w == InfNaN::neginf()) ? 0. : exp(*w - log_u))
		: ((*w == 0.) ? 0. : exp(log(*w) - log_u));

	double q_r = grat_r(b, z, log_r, eps); // = q/r of former grat1(b,z, r, &p, &q)
	double v = 0.25 / (nu * nu);
	double t2 = lnx * 0.25 * lnx;
	double j = q_r;
	double sum = j;
	double t = 1.0;
	double cn = 1.0;
	double n2 = 0.0;
	for (int n = 1; n <= n_terms_bgrat; ++n) {
		double bp2n = b + n2;
		j = (bp2n * (bp2n + 1.) * j + (z + bp2n + 1.) * t) * v;
		n2 += 2.;
		t *= t2;
		cn /= n2 * (n2 + 1.);
		int nm1 = n - 1;
		c[nm1] = cn;
		double s = 0.;
		if (n > 1) {
			double coef = b - n;
			for (int i = 1; i <= nm1; ++i) {
				s += coef * c[i - 1] * d[nm1 - i];
				coef += b;
			}
		}
		d[nm1] = bm1 * cn + s / n;
		double dj = d[nm1] * j;
		sum += dj;
		if (sum <= 0.) {
			*ierr = 3;
			return;
		}
		if (fabs(dj) <= eps * (sum + l)) {
			*ierr = 0;
			break;
		}
		else if (n == n_terms_bgrat) { // never? ; please notify R-core if seen:
			*ierr = 4;
			std::cout << "Warning: bgrat(a = " << a << ", b = " << b << ", x = " << x
				<< ") *no* convergence: NOTIFY R-core! dj = " << dj << ", rel.err = "
				<< fabs(dj) / (sum + 1) << std::endl;
		}
	} // for(n .. n_terms..)

/*                    ADD THE RESULTS TO W */

	if (log_w) // *w is in log space already:
		*w = Gamma::logspaceAdd(*w, log_u + log(sum));
	else
		*w += (u_0 ? exp(log_u + log(sum)) : u * sum);
	return;
} /* bgrat */

double Toms::brcmp1(int mu, double a, double b, double x, double y, bool give_log)
{
	const double const__ = 0.398942280401433; /* == 1/sqrt(2*pi); */
	double c = 0.0, t = 0.0, u = 0.0;
	double v = 0.0, z = 0.0, b0 = 0.0, apb = 0.0;

	double a0 = std::min(a, b);
	if (a0 < 8.) {
		double lnx, lny;
		if (x <= .375) {
			lnx = log(x);
			lny = alnrel(-x);
		}
		else if (y > .375) {
			// L11:
			lnx = log(x);
			lny = log(y);
		}
		else {
			lnx = alnrel(-y);
			lny = log(y);
		}

		// L20:
		z = a * lnx + b * lny;
		if (a0 >= 1.) {
			z -= betaln(a, b);
			return esum(mu, z, give_log);
		}
		// else :
		/* ----------------------------------------------------------------------- */
		/*              PROCEDURE FOR A < 1 OR B < 1 */
		/* ----------------------------------------------------------------------- */
		// L30:
		b0 = std::max(a, b);
		if (b0 >= 8.) {
			/* L80:                  ALGORITHM FOR b0 >= 8 */
			u = gamln1(a0) + algdiv(a0, b0);
			return give_log
				? log(a0) + esum(mu, z - u, true)
				: a0 * esum(mu, z - u, false);

		}
		else if (b0 <= 1.) {
			//                   a0 < 1, b0 <= 1
			double ans = esum(mu, z, give_log);
			if (ans == (give_log ? InfNaN::neginf() : 0.))
				return ans;

			apb = a + b;
			if (apb > 1.) {
				// L40:
				u = a + b - 1.;
				z = (gam1(u) + 1.) / apb;
			}
			else {
				z = gam1(apb) + 1.;
			}
			// L50:
			c = give_log
				? log1p(gam1(a)) + log1p(gam1(b)) - log(z)
				: (gam1(a) + 1.) * (gam1(b) + 1.) / z;
			return give_log
				? ans + log(a0) + c - log1p(a0 / b0)
				: ans * (a0 * c) / (a0 / b0 + 1.);
		}

		u = gamln1(a0);
		int n = (int)(b0 - 1.);
		if (n >= 1) {
			c = 1.;
			for (int i = 1; i <= n; ++i) {
				b0 += -1.;
				c *= b0 / (a0 + b0);
				/* L61: */
			}
			u += log(c); // TODO?: log(c) = log( prod(...) ) =  sum( log(...) )
		}
		// L70:
		z -= u;
		b0 += -1.;
		apb = a0 + b0;
		if (apb > 1.) {
			// L71:
			t = (gam1(apb - 1.) + 1.) / apb;
		}
		else {
			t = gam1(apb) + 1.;
		}

		// L72:
		return give_log
			? log(a0) + esum(mu, z, true) + log1p(gam1(b0)) - log(t) // TODO? log(t) = log1p(..)
			: a0 * esum(mu, z, false) * (gam1(b0) + 1.) / t;

	}
	else {
		double h, x0, y0, lambda;
		if (a > b) {
			// L101:
			h = b / a;
			x0 = 1. / (h + 1.);// => lx0 := log(x0) = 0 - log1p(h)
			y0 = h / (h + 1.);
			lambda = (a + b) * y - b;
		}
		else {
			h = a / b;
			x0 = h / (h + 1.);  // => lx0 := log(x0) = - log1p(1/h)
			y0 = 1. / (h + 1.);
			lambda = a - (a + b) * x;
		}
		double lx0 = -log1p(b / a); // in both cases

		double e = -lambda / a;
		if (fabs(e) > 0.6) {
			// L111:
			u = e - log(x / x0);
		}
		else {
			u = rlog1(e);
		}

		e = lambda / b;
		if (fabs(e) > 0.6) {
			v = e - log(y / y0);
		}
		else {
			v = rlog1(e);
		}

		z = esum(mu, -(a * u + b * v), give_log);
		return give_log
			? log(const__) + (log(b) + lx0) / 2. + z - bcorr(a, b)
			: const__ * sqrt(b * x0) * z * exp(-bcorr(a, b));
	}
} /* brcmp1 */

double Toms::brcomp(double a, double b, double x, double y, bool log_p)
{
	const double const__ = 0.398942280401433; /* == 1/sqrt(2*pi); */
	/* R has  M_1_SQRT_2PI , and M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938.. */
	int n = 0;
	double c = 0.0, e = 0.0, u = 0.0;
	double v = 0.0, z = 0.0, a0 = 0.0, b0 = 0.0, apb = 0.0;

	if (x == 0. || y == 0.) {
		return log_p ? InfNaN::neginf() : 0.0;
	}
	a0 = std::min(a, b);
	if (a0 < 8.) {
		double lnx, lny;
		if (x <= 0.375) {
			lnx = log(x);
			lny = alnrel(-x);
		}
		else {
			if (y > .375) {
				lnx = log(x);
				lny = log(y);
			}
			else {
				lnx = alnrel(-y);
				lny = log(y);
			}
		}

		z = a * lnx + b * lny;
		if (a0 >= 1.) {
			z -= betaln(a, b);
			return log_p ? z : exp(z);	// R_D_exp(z);
		}

		b0 = std::max(a, b);
		if (b0 >= 8.) { /* L80: */
			u = gamln1(a0) + algdiv(a0, b0);
			return (log_p ? log(a0) + (z - u) : a0 * exp(z - u));
		}

		if (b0 <= 1.) { /*		algorithm for max(a,b) = b0 <= 1 */

			double e_z = log_p ? z : exp(z);

			if (!log_p && e_z == 0.) /* exp() underflow */
				return 0.;

			apb = a + b;
			if (apb > 1.) {
				u = a + b - 1.;
				z = (gam1(u) + 1.) / apb;
			}
			else {
				z = gam1(apb) + 1.;
			}

			c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;
			/* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
			return (log_p
				? e_z + log(a0 * c) - log1p(a0 / b0)
				: e_z * (a0 * c) / (a0 / b0 + 1.));
		}

		u = gamln1(a0);
		n = (int)(b0 - 1.);
		if (n >= 1) {
			c = 1.;
			for (int i = 1; i <= n; ++i) {
				b0 += -1.;
				c *= b0 / (a0 + b0);
			}
			u = log(c) + u;
		}
		z -= u;
		b0 += -1.;
		apb = a0 + b0;
		double t;
		if (apb > 1.) {
			u = a0 + b0 - 1.;
			t = (gam1(u) + 1.) / apb;
		}
		else {
			t = gam1(apb) + 1.;
		}

		return (log_p
			? log(a0) + z + log1p(gam1(b0)) - log(t)
			: a0 * exp(z) * (gam1(b0) + 1.) / t);

	}
	else {
		double h = 0.0, x0 = 0.0, y0 = 0.0, lambda = 0.0;
		if (a <= b) {
			h = a / b;
			x0 = h / (h + 1.);
			y0 = 1. / (h + 1.);
			lambda = a - (a + b) * x;
		}
		else {
			h = b / a;
			x0 = 1. / (h + 1.);
			y0 = h / (h + 1.);
			lambda = (a + b) * y - b;
		}

		e = -lambda / a;
		if (fabs(e) > .6)
			u = e - log(x / x0);
		else
			u = rlog1(e);

		e = lambda / b;
		if (fabs(e) <= .6)
			v = rlog1(e);
		else
			v = e - log(y / y0);

		z = log_p ? -(a * u + b * v) : exp(-(a * u + b * v));

		return log_p
			? -M_LN_SQRT_2PI + .5 * log(b * x0) + z - bcorr(a, b)
			: const__ * sqrt(b * x0) * z * exp(-bcorr(a, b));
	}
} /* brcomp */

double Toms::bfrac(double a, double b, double x, double y, double lambda,
	double eps, bool log_p)
{
	double e = 0.0, t = 0.0, w = 0.0, r0 = 0.0;
	double beta = 0.0, alpha = 0.0, brc = 0.0;

	if (!std::isfinite(lambda)) return InfNaN::nan();// TODO: can return 0 or 1 (?)

	brc = brcomp(a, b, x, y, log_p);
	if (std::isnan(brc)) { // e.g. from   L <- 1e308; pnbinom(L, L, mu = 5)
		InfNaN::nan(); // TODO: could we know better?
	}
	if (!log_p && brc == 0.0) {
		return 0.0;
	}

	double c = lambda + 1.0;
	double c0 = b / a;
	double c1 = 1.0 / a + 1.0;
	double yp1 = y + 1.0;

	double n = 0.0;
	double p = 1.0;
	double s = a + 1.0;
	double an = 0.0;
	double bn = 1.0;
	double anp1 = 1.0;
	double bnp1 = c / c1;
	double r = c1 / c;

	/*        CONTINUED FRACTION CALCULATION */

	do {
		n += 1.;
		t = n / a;
		w = n * (b - n) * x;
		e = a / s;
		alpha = p * (p + c0) * e * e * (w * x);
		e = (t + 1.) / (c1 + t + t);
		beta = n + w / s + e * (c + n * yp1);
		p = t + 1.;
		s += 2.;

		/* update an, bn, anp1, and bnp1 */

		t = alpha * an + beta * anp1;	an = anp1;	anp1 = t;
		t = alpha * bn + beta * bnp1;	bn = bnp1;	bnp1 = t;

		r0 = r;
		r = anp1 / bnp1;

		if (fabs(r - r0) <= eps * r)
			break;

		/* rescale an, bn, anp1, and bnp1 */

		an /= bnp1;
		bn /= bnp1;
		anp1 = r;
		bnp1 = 1.;
	} while (n < 10000);// arbitrary; had '1' --> infinite loop for  lambda = Inf

	if (n >= 10000 && fabs(r - r0) > eps * r)
		std::cout << "Warning: bfrac(a = " << a << ", b = " << b << " x = " << x
			<< ", y = " << y << ", lambda = " << lambda << ") did *not* converge"
			<< " (in 10000 steps)" << std::endl;
	return (log_p ? brc + log(r) : brc * r);
} /* bfrac */

double Toms::bup(double a, double b, double x, double y, int n, double eps, bool give_log)
{
	double ret_val = 0.0;
	int k = 0, mu = 0;
	double d = 0.0, l = 0.0;


	double apb = a + b;
	double ap1 = a + 1.0;
	if (n > 1 && a >= 1. && apb >= ap1 * 1.1) {
		mu = (int)fabs(exparg(1));
		k = (int)exparg(0);
		if (mu > k)
			mu = k;
		d = exp(-(double)mu);
	}
	else {
		mu = 0;
		d = 1.0;
	}

	/* L10: */
	ret_val = give_log
		? brcmp1(mu, a, b, x, y, true) - log(a)
		: brcmp1(mu, a, b, x, y, false) / a;
	if (n == 1 ||
		(give_log && ret_val == InfNaN::neginf()) || (!give_log && ret_val == 0.))
		return ret_val;

	int nm1 = n - 1;
	double w = d;

	k = 0;
	if (b > 1.) {
		if (y > 1e-4) {
			double r = (b - 1.) * x / y - a;
			if (r >= 1.)
				k = (r < nm1) ? (int)r : nm1;
		}
		else
			k = nm1;

		for (int i = 0; i < k; ++i) {
			l = (double)i;
			d *= (apb + l) / (ap1 + l) * x;
			w += d;
		}
	}

	for (int i = k; i < nm1; ++i) {
		l = (double)i;
		d *= (apb + l) / (ap1 + l) * x;
		w += d;
		if (d <= eps * w) /* relativ convergence (eps) */
			break;
	}

	if (give_log) {
		ret_val += log(w);
	}
	else
		ret_val *= w;
	return ret_val;
} /* bup */

double Toms::bpser(double a, double b, double x, double eps, bool log_p)
{
	int m = 0;
	double ans = 0.0, c = 0.0, t = 0.0;
	double u = 0.0, z = 0.0, b0 = 0.0, apb = 0.0;

	if (x == 0.0) {
		return log_p ? InfNaN::neginf() : 0.0;
	}

	double a0 = std::min(a, b);
	if (a0 >= 1.) { /*		 ------	 1 <= a0 <= b0  ------ */
		z = a * log(x) - betaln(a, b);
		ans = log_p ? z - log(a) : exp(z) / a;
	}
	else {
		b0 = std::max(a, b);

		if (b0 < 8.) {

			if (b0 <= 1.) { /*	 ------	 a0 < 1	 and  b0 <= 1  ------ */

				if (log_p) {
					ans = a * log(x);
				}
				else {
					ans = pow(x, a);
					if (ans == 0.) /* once underflow, always underflow .. */
						return ans;
				}
				apb = a + b;
				if (apb > 1.) {
					u = a + b - 1.;
					z = (gam1(u) + 1.) / apb;
				}
				else {
					z = gam1(apb) + 1.;
				}
				c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;

				if (log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
					ans += log(c * (b / apb));
				else
					ans *= c * (b / apb);

			}
			else { /* 	------	a0 < 1 < b0 < 8	 ------ */

				u = gamln1(a0);
				m = (int)(b0 - 1.);
				if (m >= 1) {
					c = 1.;
					for (int i = 1; i <= m; ++i) {
						b0 += -1.;
						c *= b0 / (a0 + b0);
					}
					u += log(c);
				}

				z = a * log(x) - u;
				b0 += -1.; // => b0 in (0, 7)
				apb = a0 + b0;
				if (apb > 1.) {
					u = a0 + b0 - 1.;
					t = (gam1(u) + 1.) / apb;
				}
				else {
					t = gam1(apb) + 1.;
				}

				if (log_p) /* FIXME? potential for improving log(t) */
					ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t);
				else
					ans = exp(z) * (a0 / a) * (gam1(b0) + 1.) / t;
			}
		}
		else { /* 		------  a0 < 1 < 8 <= b0  ------ */

			u = gamln1(a0) + algdiv(a0, b0);
			z = a * log(x) - u;

			if (log_p)
				ans = z + log(a0 / a);
			else
				ans = a0 / a * exp(z);
		}
	}
	if ((ans == (log_p ? InfNaN::neginf() : 0.0)) || (!log_p && a <= eps * 0.1)) {
		return ans;
	}

	double tol = eps / a;
	double n = 0.0;
	double sum = 0.0;
	double w = 0.0;
	c = 1.0;
	do { // sum is alternating as long as n < b (<==> 1 - b/n < 0)
		n += 1.;
		c *= (0.5 - b / n + 0.5) * x;
		w = c / (a + n);
		sum += w;
	} while (n < 1e7 && fabs(w) > tol);

	if (fabs(w) > tol) {
		if ((log_p && !(a * sum > -1. && fabs(log1p(a * sum)) < eps * fabs(ans))) ||
			(!log_p && fabs(a * sum + 1.) != 1.))
				std::cout << "Warngin: bpser(a = " << a << ", b = " << b << ", x = "
					<< x << ",...) did not converge (n=1e7, |w|/tol = " << fabs(w) / tol
					<< " > 1; A = " << ans << std::endl;
	}

	if (log_p) {
		if (a * sum > -1.0) ans += log1p(a * sum);
		else {
			if (ans > InfNaN::neginf())
				std::cout << "Warning: pbeta(*, log.p=TRUE) -> bpser(a = " << a << ", b = "
					<< b << ", x = " << x << ",...) underflow to -Inf" << std::endl;
			ans = InfNaN::neginf();
		}
	}
	else if (a * sum > -1.)
		ans *= (a * sum + 1.);
	else // underflow to
		ans = 0.;
	return ans;
} /* bpser */

double Toms::apser(double a, double b, double x, double eps)
{
	const double g = 0.577215664901533;

	double c = 0.0, aj = 0.0;
	double bx = b * x;

	double t = x - bx;
	if (b * eps <= 0.02)
		c = log(x) + psi(b) + g + t;
	else // b > 2e13 : psi(b) ~= log(b)
		c = log(bx) + g + t;

	double tol = eps * 5. * fabs(c);
	double j = 1.0;
	double s = 0.;
	do {
		j += 1.;
		t *= x - bx / j;
		aj = t / j;
		s += aj;
	} while (fabs(aj) > tol);

	return -a * (c + s);
} /* apser */

double Toms::fpser(double a, double b, double x, double eps, bool log_p)
{
	double ans = 0.0, c = 0.0, t = 0.0;

	/* SET  ans := x^a : */
	if (log_p) {
		ans = a * log(x);
	}
	else if (a > eps * 0.001) {
		t = a * log(x);
		if (t < exparg(1)) { /* exp(t) would underflow */
			return 0.;
		}
		ans = exp(t);
	}
	else
		ans = 1.;

	if (log_p)
		ans += log(b) - log(a);
	else
		ans *= b / a;

	double tol = eps / a;
	double an = a + 1.0;
	t = x;
	double s = t / an;
	do {
		an += 1.;
		t = x * t;
		c = t / an;
		s += c;
	} while (fabs(c) > tol);

	if (log_p)
		ans += log1p(a * s);
	else
		ans *= a * s + 1.;
	return ans;
} /* fpser */