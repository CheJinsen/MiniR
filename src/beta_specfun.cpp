#include "specfun.h"

using namespace SpecialFunctions;


double Beta::lbeta(const double a, const double b)
{
	double corr = 0.0, p = 0.0, q = 0.0;

	if (std::isnan(a) || std::isnan(b))
		return a + b;

	p = q = a;
	if (b < p) p = b;/* := min(a,b) */
	if (b > q) q = b;/* := max(a,b) */

	/* both arguments must be >= 0 */
	if (p < 0) {
		return InfNaN::nan();
	}
	else if (p == 0) {
		return InfNaN::posinf();
	}
	else if (!std::isfinite(q)) { /* q == +Inf */
		return InfNaN::neginf();
	}

	if (p >= 10) {
		/* p and q are big. */
		corr = Gamma::lgammacor(p) + Gamma::lgammacor(q) - Gamma::lgammacor(p + q);
		return log(q) * -0.5 + M_LN_SQRT_2PI + corr +
			(p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
	}
	else if (q >= 10) {
		/* p is small, but q is big. */
		corr = Gamma::lgammacor(q) - Gamma::lgammacor(p + q);
		return Gamma::lgammafn(p) + corr + p - p * log(p + q) +
			(q - 0.5) * log1p(-p / (p + q));
	}
	else {
		/* p and q are small: p <= q < 10. */
		/* R change for very small args */
		if (p < 1e-306)
			return lgamma(p) + (lgamma(q) - lgamma(p + q));
		else
			return log(Gamma::gammafn(p) *
				(Gamma::gammafn(q) / Gamma::gammafn(p + q)));
	}
}