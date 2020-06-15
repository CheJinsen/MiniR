
#pragma once

#include "base.h"

double dnorm(double x, double mean = 0.0,
		double sigma = 1.0, bool log_p = false);

double pnorm(double x, double mean = 0.0,
		double sigma = 1.0, bool lower_tail = true, bool log_p = false);

void pnorm_both(double x, double& cum, double& ccum, bool i_tail, bool log_p);
