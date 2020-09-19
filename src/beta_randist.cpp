/*
 * This file is part of MiniR.
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 2000-2018 The R Core Team
 * Copyright (C) 2001-2016 The R Foundation
 * Copyright (C) 2020 Jinsen Che
 *
 * MiniR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MiniR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MiniR. If not, see <https://www.gnu.org/licenses/>.
 */

#include "randist.h"

using namespace Randist;

double Beta::cdfRaw(double x, double a, double b,
	bool lower_tail, bool log_p)
{
	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	if (a == 0 || b == 0 || !std::isfinite(a) || !std::isfinite(b)) {
		if (a == 0 && b == 0)
			return (log_p ? -M_LN2 : 0.5);
		if (a == 0 || a / b == 0)
			return R_DT_1;
		if (b == 0 || b / a == 0)
			return R_DT_0;
		if (x < 0.5) return R_DT_0; else return R_DT_1;
	}

	double x1 = 0.5 - x + 0.5, w, wc;
	int ierr;

	Toms::bratio(a, b, x, x1, &w, &wc, &ierr, log_p); // acmtoms.cpp

	if (ierr && ierr != 11 && ierr != 14)
		std::cout << "Warning: pbeta_raw(" << x << ", a = " << a << ", b = "
		<< b << ", ..) -> bratio() gave error code " << ierr << std::endl;
	return lower_tail ? w : wc;
}

double Beta::cdf(double x, double a, double b, bool lower_tail, bool log_p)
{
	if (std::isnan(x) || std::isnan(a) || std::isnan(b)) return x + a + b;

	if (a < 0 || b < 0) return InfNaN::nan();

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;
	if (x <= 0)
		return R_DT_0;
	if (x >= 1)
		return R_DT_1;

	return cdfRaw(x, a, b, lower_tail, log_p);
}

double Beta::pdf(double x, double a, double b, bool give_log)
{
	if (std::isnan(x) || std::isnan(a) || std::isnan(b)) return x + a + b;
	double R_D__0 = give_log ? InfNaN::neginf() : 0.0;

	if (a < 0 || b < 0) return InfNaN::nan();
	if (x < 0 || x > 1) return R_D__0;

	if (a == 0 || b == 0 || !std::isfinite(a) || !std::isfinite(b)) {
		if (a == 0 && b == 0) {
			if (x == 0 || x == 1)
				return InfNaN::posinf();
			else
				return R_D__0 ;
		}
		if (a == 0 || a / b == 0) {
			if (x == 0)
				return InfNaN::posinf();
			else
				return R_D__0;
		}
		if (b == 0 || b / a == 0) {
			if (x == 1)
				return InfNaN::posinf();
			else
				return R_D__0;
		}
		if (x == 0.5)
			return InfNaN::posinf();
		else
			return R_D__0;
	}

	if (x == 0) {
		if (a > 1) return R_D__0;
		if (a < 1) return InfNaN::posinf();
		return give_log ? log(b) : (b);
	}
	if (x == 1) {
		if (b > 1) return R_D__0;
		if (b < 1) return InfNaN::posinf();
		return give_log ? log(a) : (a);
	}

	double lval = 0.0;
	if (a <= 2 || b <= 2)
		lval = (a - 1) * log(x) + (b - 1) * log1p(-x) - SpecialFunctions::Beta::lbeta(a, b);
	else
		lval = log(a + b - 1) + Base::binomialPdfRaw(a - 1, a + b - 2, x, 1 - x, true);

	return give_log ? lval : exp(lval);
}

double Beta::quantile(double alpha, double p, double q,
	bool lower_tail, bool log_p)
{
	if (std::isnan(p) || std::isnan(q) || std::isnan(alpha))
		return p + q + alpha;

	if (p < 0. || q < 0.) return InfNaN::nan();

	double qbet[2];// = { qbeta(), 1 - qbeta() }
	quantileRaw(alpha, p, q, lower_tail, log_p, -1, -5, 4, qbet);
	return qbet[0];
}

double Beta::powDi(double x, int n)
{
	double vpow = 1.0;

	if (std::isnan(x)) return x;
	if (n != 0) {
		if (!std::isfinite(x)) return pow(x, (double)n);
		if (n < 0) { n = -n; x = 1 / x; }
		for (;;) {
			if (n & 01) vpow *= x;
			if (n >>= 1) x *= x; else break;
		}
	}
	return vpow;
}

void Beta::quantileRaw(double alpha, double p, double q,
	bool lower_tail, bool log_p, int swap_01, double log_q_cut,
	int n_N, double* qb)
{
	bool swap_choose = (swap_01 == -1);
	bool swap_tail = false;
	bool log_ = false;
	bool give_log_q = (log_q_cut == InfNaN::posinf());
	bool use_log_x = give_log_q; // or u < log_q_cut  below
	bool warned = false;
	bool add_N_step = true;

	double R_D__0 = log_p ? InfNaN::neginf() : 0.0;
	double R_D__1 = log_p ? 0.0 : 1.0;
	double R_DT_0 = lower_tail ? R_D__0 : R_D__1;
	double R_DT_1 = lower_tail ? R_D__1 : R_D__0;

	const double DBL_very_MIN = DBL_MIN / 4.0;
	const double DBL_log_v_MIN = M_LN2 * (DBL_MIN_EXP - 2.0);
	constexpr double p_lo = 3e-308;
	constexpr double DBL_1__eps = 0x1.fffffffffffffp-1;
	double wprev = 0.0, prev = 1.0, adj = 1.0;
	double g = 0.0;

	if (alpha == R_DT_0) {
#define return_q_0						\
	if(give_log_q) { qb[0] = InfNaN::neginf(); qb[1] = 0; }	\
	else {           qb[0] = 0;         qb[1] = 1; }	\
	return

		return_q_0;
	}
	if (alpha == R_DT_1) {
#define return_q_1						\
	if(give_log_q) { qb[0] = 0; qb[1] = InfNaN::neginf(); }	\
	else {           qb[0] = 1; qb[1] = 0;         }	\
	return

		return_q_1;
	}

	if ((log_p && alpha > 0) ||
		(!log_p && (alpha < 0 || alpha > 1))) {
		std::cout << "Argument out of domain in beta_randist.cpp" << std::endl;
		qb[0] = qb[1] = InfNaN::nan();
		return;
	}

	double R_D_half = log_p ? -M_LN2 : 0.5;
	if (p == 0 || q == 0 || !std::isfinite(p) || !std::isfinite(q)) {
		if (p == 0 && q == 0) {
			if (alpha < R_D_half) { return_q_0; }
			if (alpha > R_D_half) { return_q_1; }
			
#define return_q_half					\
	    if(give_log_q) qb[0] = qb[1] = -M_LN2;	\
	    else	   qb[0] = qb[1] = 0.5;		\
	    return

			return_q_half;
		}
		else if (p == 0 || p / q == 0) {
			return_q_0;
		}
		else if (q == 0 || q / p == 0) {
			return_q_1;
		}
		return_q_half;
	}

	double p_ = log_p ? (lower_tail ? exp(alpha) : -expm1(alpha))
		: (lower_tail ? alpha : (0.5 - alpha + 0.5));
	double logbeta = SpecialFunctions::Beta::lbeta(p, q);

	swap_tail = (swap_choose) ? (p_ > 0.5) : swap_01;
	double a = 0.0, la = 0.0, pp = 0.0;
	double qq = 0.0;
	if (swap_tail) {
		a = log_p ? (lower_tail ? -expm1(alpha) : exp(alpha))
			: (lower_tail ? (0.5 - alpha + 0.5) : alpha);

		double t1 = log_p ? alpha : log(alpha);
		double t2 = alpha > -M_LN2 ? log(-expm1(alpha)) : log1p(-exp(alpha));
		double t3 = log_p ? t2 : log1p(-alpha);
		la = lower_tail ? t3 : t1;
		pp = q; qq = p;
	}
	else {
		double t0 = log_p ? alpha : log(alpha);
		double t1 = alpha > -M_LN2 ? log(-expm1(alpha)) : log1p(-exp(alpha));
		double t2 = log_p ? t1 : log1p(-alpha);
		
		a = p_;
		la = lower_tail ? t0 : t2;
		pp = p; qq = q;
	}

	constexpr auto acu_min = 1e-300;
	constexpr auto fpu = 3e-308;
	constexpr auto p_hi = 1-2.22e-16;
	constexpr auto const1 = 2.30753;
	constexpr auto const2 = 0.27061;
	constexpr auto const3 = 0.99229;
	constexpr auto const4 = 0.04481;

	double acu = std::max(acu_min, pow(10., -13. - 2.5 / (pp * pp) - 0.5 / (a * a)));

	double tx = 0.0;
	double u0 = (la + log(pp) + logbeta) / pp;
	const double log_eps_c = M_LN2 * (1. - DBL_MANT_DIG);
	double r = pp * (1. - qq) / (pp + 1.);
	double t = 0.2;
	double u = 0.0;
	double xinbta = 0.0;
	double u_n = 0.0;
	double y = 0.0;

	bool bad_u = false;
	bool bad_init = false;

	double s = 0.0;
	double h = 0.0;
	double w = 0.0;

	if (M_LN2 * DBL_MIN_EXP < u0 && u0 < -0.01 &&
		u0 < (t * log_eps_c - log(fabs(pp * (1.0 - qq) * (2.0 - qq) /
		(2.0 * (pp + 2.0))))) / 2.0)
	{
		r = r * exp(u0);// = r*x0
		if (r > -1.) {
			u = u0 - log1p(r) / pp;
		}
		else {
			u = u0;
		}
		tx = xinbta = exp(u);
		use_log_x = true;
		goto L_Newton;
	}

	r = sqrt(-2 * la);
	y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);

	if (pp > 1 && qq > 1) {
		r = (y * y - 3.) / 6.;
		s = 1. / (pp + pp - 1.);
		t = 1. / (qq + qq - 1.);
		h = 2. / (s + t);
		w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));

		if (w > 300) {
			t = w + w + log(qq) - log(pp);
			u = (t <= 18) ? -log1p(exp(t)) : -t - exp(-t);
			xinbta = exp(u);
		}
		else {
			xinbta = pp / (pp + qq * exp(w + w));
			u = -log1p(qq / pp * exp(w + w));
		}
	}
	else {
		r = qq + qq;
		t = 1. / (3. * sqrt(qq));
		t = r * powDi(1. + t * (-t + y), 3);
		s = 4. * pp + r - 2.0;

		if (t == 0 || (t < 0. && s >= t)) {
			double l1ma = 0.0;
			double t1 = alpha > -M_LN2 ? log(-expm1(alpha)) : log1p(-exp(alpha));
			double t2 = log_p ? t1 : log1p(-alpha);
			double t3 = log_p ? alpha : log(alpha);
			if (swap_tail)
				l1ma = lower_tail ? t3 : t2;
			else
				l1ma = lower_tail ? t2 : t3;

			double xx = (l1ma + log(qq) + logbeta) / qq;
			if (xx <= 0.) {
				xinbta = -expm1(xx);
				u = xx > -M_LN2 ? log(-expm1(xx)) : log1p(-exp(xx));
			}
			else {
				xinbta = 0; u = InfNaN::nan();
			}
		}
		else {
			t = s / t;
			if (t <= 1.) {
				u = (la + log(pp) + logbeta) / pp;
				xinbta = exp(u);
			}
			else {
				xinbta = 1. - 2. / (t + 1.);
				u = log1p(-2. / (t + 1.));
			}
		}
	}

	if (swap_choose &&
		((swap_tail && u >= -exp(log_q_cut)) || // ==> "swap back"
			(!swap_tail && u >= -exp(4 * log_q_cut) && pp / qq < 1000.))) {

		double t1 = alpha > -M_LN2 ? log(-expm1(alpha)) : log1p(-exp(alpha));
		double t2 = log_p ? t1 : log1p(-alpha);
		double t3 = log_p ? alpha : log(alpha);
		double t4 = lower_tail ? (0.5 - alpha + 0.5) : alpha;

		swap_tail = !swap_tail;
		if (swap_tail) {
			a = log_p ? (lower_tail ? -expm1(alpha) : exp(alpha)) : t4;
			la = lower_tail ? t2 : t3;
			pp = q; qq = p;
		}
		else {
			a = p_;
			la = lower_tail ? t3 : t2;
			pp = p; qq = q;
		}
		u = u > -M_LN2 ? log(-expm1(u)) : log1p(-exp(u));
		xinbta = exp(u);
	}

	if (!use_log_x)
		use_log_x = (u < log_q_cut);
	
	bad_u = !std::isfinite(u);
	bad_init = (bad_u || xinbta > p_hi);

	u_n = 1.0; // -Wall
	tx = xinbta; // keeping "original initial x" (for now)

	if (bad_u || u < log_q_cut) {
		w = cdfRaw(DBL_very_MIN, pp, qq, true, log_p);
		if (w > (log_p ? la : a)) {
			if (log_p || fabs(w - a) < fabs(0 - a)) {
				tx = DBL_very_MIN;
				u_n = DBL_log_v_MIN;
			}
			else {
				tx = 0.;
				u_n = InfNaN::neginf();
			}
			use_log_x = log_p; add_N_step = false;
			goto L_return;
		}
		else {
			if (u < DBL_log_v_MIN) {
				u = DBL_log_v_MIN;
				xinbta = DBL_very_MIN;
			}
		}
	}

	if (bad_init && !(use_log_x && tx > 0)) {
		if (u == InfNaN::neginf()) {
			u = M_LN2 * DBL_MIN_EXP;
			xinbta = DBL_MIN;
		}
		else {
			xinbta = (xinbta > 1.1) ? 0.5 : ((xinbta < p_lo) ? exp(u) : p_hi);
			if (bad_u)
				u = log(xinbta);
		}
	}

L_Newton:

	r = 1 - pp;
	t = 1 - qq;

	if (use_log_x) {
		for (int i_pb = 0; i_pb < 1000; i_pb++) {
			y = cdfRaw(xinbta, pp, qq, true, true);
			double temp = u > -M_LN2 ? log(-expm1(u)) : log1p(-exp(u));
			w = (y == InfNaN::neginf())
				? 0. : (y - la) * exp(y - u + logbeta + r * u + t * temp);
			if (!std::isfinite(w))
				break;
			if (i_pb >= n_N && w * wprev <= 0.)
				prev = std::max(fabs(adj), fpu);
			g = 1;
			for (int i_inn = 0; i_inn < 1000; i_inn++) {
				adj = g * w;
				if (fabs(adj) < prev) {
					u_n = u - adj;
					if (u_n <= 0.) {
						if (prev <= acu || fabs(w) <= acu) {
							goto L_converged;
						}
						break;
					}
				}
				g /= 3;
			}
			
			double D = std::min(fabs(adj), fabs(u_n - u));
			if (D <= 4e-16 * fabs(u_n + u))
				goto L_converged;
			u = u_n;
			xinbta = exp(u);
			wprev = w;
		} 
	}
	else {
		for (int i_pb = 0; i_pb < 1000; i_pb++) {
			y = cdfRaw(xinbta, pp, qq, true, log_p);

			if (!std::isfinite(y) && !(log_p && y == InfNaN::neginf())) {
				std::cout << "Argument out of domain in beta_randist.cpp" << std::endl;
				qb[0] = qb[1] = InfNaN::nan();
				return;
			}

			w = log_p
				? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
				: (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
			if (i_pb >= n_N && w * wprev <= 0.)
				prev = std::max(fabs(adj), fpu);
			g = 1.0;
			for (int i_inn = 0; i_inn < 1000; i_inn++) {
				adj = g * w;
				if (i_pb < n_N || fabs(adj) < prev) {
					tx = xinbta - adj;
					if (0. <= tx && tx <= 1.) {
						if (prev <= acu || fabs(w) <= acu) {
							goto L_converged;
						}
						if (tx != 0. && tx != 1)
							break;
					}
				}
				g /= 3;
			}
			if (fabs(tx - xinbta) <= 4e-16 * (tx + xinbta))
				goto L_converged;
			xinbta = tx;
			if (tx == 0)
				break;
			wprev = w;
		}

	}
	warned = true;
	std::cout << "Full precision may not have been achieved in beta quantile." << std::endl;

L_converged:
	log_ = log_p || use_log_x;

	if ((log_ && y == InfNaN::neginf()) || (!log_ && y == 0)) {
		w = cdfRaw(DBL_very_MIN, pp, qq, true, log_);
		if (log_ || fabs(w - a) <= fabs(y - a)) {
			tx = DBL_very_MIN;
			u_n = DBL_log_v_MIN;
		}
		add_N_step = false;
	}
	else if (!warned && (log_ ? fabs(y - la) > 3 : fabs(y - a) > 1e-4)) {
		if (!(log_ && y == InfNaN::neginf() &&
			cdfRaw(DBL_1__eps, pp, qq, true, true) > la + 2))
			//MATHLIB_WARNING2( // low accuracy for more platform independent output:
			//	"qbeta(a, *) =: x0 with |pbeta(x0,*%s) - alpha| = %.5g is not accurate",
			//	(log_ ? ", log_" : ""), fabs(y - (log_ ? la : a)));
			std::cout << "Warning: beta quantile function is not accurate." << std::endl;
	}
L_return:
	if (give_log_q) {
		if (!use_log_x)
			std::cout << "Warning: Beta::quantile() L_return, u_n = " << u_n
				<< "; give_log_q = TRUE but use_log_x = FALSE -- please report!" << std::endl;
		double r = u_n > -M_LN2 ? log(-expm1(u_n)) : log1p(-exp(u_n));
		if (swap_tail) {
			qb[0] = r;	 qb[1] = u_n;
		}
		else {
			qb[0] = u_n; qb[1] = r;
		}
	}
	else {
		if (use_log_x) {
			if (add_N_step) {
				xinbta = exp(u_n);
				y = cdfRaw(xinbta, pp, qq, true, log_p);
				w = log_p
					? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
					: (y - a) * exp(logbeta + r * log(xinbta) + t * log1p(-xinbta));
				tx = xinbta - w;
			}
			else {
				if (swap_tail) {
					qb[0] = -expm1(u_n); qb[1] = exp(u_n);
				}
				else {
					qb[0] = exp(u_n); qb[1] = -expm1(u_n);
				}
				return;
			}
		}
		if (swap_tail) {
			qb[0] = 1 - tx; qb[1] = tx;
		}
		else {
			qb[0] = tx;	qb[1] = 1 - tx;
		}
	}
	return;
}

void Beta::showBet(double aa, double beta, double u1, double& v, double& w)
{
	constexpr double expmax = DBL_MAX_EXP * M_LN2;
	v = beta * log(u1 / (1.0 - u1));
	if (v <= expmax) {
		w = aa * exp(v);
		if (!std::isfinite(w)) w = DBL_MAX;
	}
	else {
		w = DBL_MAX;
	}
}

double Beta::rand(const double aa, const double bb)
{
    if (std::isnan(aa) || std::isnan(bb) || aa < 0.0 || bb < 0.0) {
		return InfNaN::nan();
    }
    
    if (!std::isfinite(aa) && !std::isfinite(bb)) {
		return 0.5;
    }
    
    if (aa == 0.0 && bb == 0.0) {
		return (Uniform::rand() < 0.5) ? 0.0 : 1.;
    }
 
    if (!std::isfinite(aa) || bb == 0.0) {
    	return 1.0;
    }
    
    if (!std::isfinite(bb) || aa == 0.0) {
    	return 0.0;
    }

    double u1 = 0.0, u2 = 0.0, v = 0.0;
    double w = 0.0, y = 0.0, z = 0.0;
    
    double beta = 0.0;
    bool qsame = (-1.0 == aa) && (-1.0 == bb);
    double a = std::min(aa, bb);
    double b = std::max(aa, bb);
    double alpha = a + b;

    if (a <= 1.0) {
    	double k1 = 0.0, k2 = 0.0;

		if (!qsame) {
		    beta = 1.0 / a;
		    double delta = 1.0 + b - a;
		    k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
		    k2 = 0.25 + (0.5 + 0.25 / delta) * a;
		}
		
		for(;;) {
		    u1 = Uniform::rand();
		    u2 = Uniform::rand();
		    if (u1 < 0.5) {
				y = u1 * u2;
				z = u1 * y;
				if (0.25 * u2 + z - y >= k1)
				    continue;
		    }
		    else {
				z = u1 * u1 * u2;
				if (z <= 0.25) {
				    showBet(b, beta, u1, v, w);
				    break;
				}
				if (z >= k2)
				    continue;
		    }

		    showBet(b, beta, u1, v, w);

		    if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
				break;
		}
		return (aa == a) ? a / (a + w) : w / (a + w);
    }
    else {
    	double gamma = 0.0;

		if (!qsame) {
		    beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
		    gamma = a + 1.0 / beta;
		}

		double r = 0.0;
		double s = 0.0;
		double t = 0.0;
		do {
		    u1 = Uniform::rand();
		    u2 = Uniform::rand();

		    showBet(a, beta, u1, v, w);

		    z = u1 * u1 * u2;
		    r = gamma * v - 1.3862944;
		    s = a + r - w;
		    if (s + 2.609438 >= 5.0 * z)
				break;
		    t = log(z);
		    if (s > t)
				break;
		} while (r + alpha * log(alpha / (b + w)) < t);

		return (aa != a) ? b / (b + w) : w / (b + w);
    }
}