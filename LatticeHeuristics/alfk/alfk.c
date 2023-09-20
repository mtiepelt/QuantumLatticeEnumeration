#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#define CDF_ERROR(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, GSL_NAN)

static double 
bisect (double x, double P, double a, double b, double xtol, double Ptol)
{
  double x0 = 0, x1 = 1, Px;

  while (fabs(x1 - x0) > xtol) {
    Px = gsl_cdf_beta_P (x, a, b);
    if (fabs(Px - P) < Ptol) {
      /* return as soon as approximation is good enough, including on
         the first iteration */
      return x;  
    } else if (Px < P) {
      x0 = x;
    } else if (Px > P) {
      x1 = x;
    }
    x = 0.5 * (x0 + x1);
  }
  return x;
}  

double
my_gsl_cdf_beta_Pinv (const double P, const double a, const double b, const unsigned int max_iterations, const double convergence_err)
{
  double x, mean;

  if (P < 0.0 || P > 1.0)
    {
      CDF_ERROR ("P must be in range 0 < P < 1", GSL_EDOM);
    }

  if (a < 0.0)
    {
      CDF_ERROR ("a < 0", GSL_EDOM);
    }

  if (b < 0.0)
    {
      CDF_ERROR ("b < 0", GSL_EDOM);
    }

  if (P == 0.0)
    {
      return 0.0;
    }

  if (P == 1.0)
    {
      return 1.0;
    }

  if (P > 0.5)
    {
      // printf("calling Qinv\n");
      return gsl_cdf_beta_Qinv (1 - P, a, b);
    }

  mean = a / (a + b);

  if (P < 0.1)
    {
      /* small x */

      double lg_ab = gsl_sf_lngamma (a + b);
      double lg_a = gsl_sf_lngamma (a);
      double lg_b = gsl_sf_lngamma (b);

      double lx = (log (a) + lg_a + lg_b - lg_ab + log (P)) / a;
      if (lx <= 0) {
        x = exp (lx);             /* first approximation */
        x *= pow (1 - x, -(b - 1) / a);   /* second approximation */
      } else {
        x = mean;
      }

      if (x > mean)
        x = mean;
    }
  else
    {
      /* Use expected value as first guess */
      x = mean;
    }

  /* Do bisection to get closer */
  x = bisect (x, P, a, b, 0.01, 0.01);

  {
    double lambda, dP, phi;
    unsigned int n = 0;

  start:
    dP = P - gsl_cdf_beta_P (x, a, b);
    phi = gsl_ran_beta_pdf (x, a, b);

    if (dP == 0.0 || n++ > max_iterations) { // Fernando: was n++ > 64
      // if (dP == 0.0) printf("issue is dP\n");
      // if (n > 128) printf("issue is n\n");
      goto end;
    }

    lambda = dP / GSL_MAX (2 * fabs (dP / x), phi);

    {
      double step0 = lambda;
      double step1 = -((a - 1) / x - (b - 1) / (1 - x)) * lambda * lambda / 2;

      double step = step0;

      if (fabs (step1) < fabs (step0))
        {
          step += step1;
        }
      else
        {
          /* scale back step to a reasonable size when too large */
          step *= 2 * fabs (step0 / step1);
        };

      if (x + step > 0 && x + step < 1)
        {
          x += step;
        }
      else
        {
          x = sqrt (x) * sqrt (mean);   /* try a new starting point */
        }

      if (fabs (step0) > convergence_err * x) // Fernando: was 1e-10
        goto start;
    }

  end:

    if (fabs(dP) > GSL_SQRT_DBL_EPSILON * P)
      {
        // printf("fabs(dP) = %g > %g = %g * %g = GSL_SQRT_DBL_EPSILON * P\n", fabs(dP), GSL_SQRT_DBL_EPSILON * P, GSL_SQRT_DBL_EPSILON, P);
        // printf("about to complain\n");
        // GSL_ERROR_VAL("inverse failed to converge", GSL_EFAILED, GSL_NAN);
        return GSL_NAN;
      }

    return x;
  }
}

// int main ()
// {
// 	double alpha = ((double)1.) / (1<<24);
// 	int k = 540;
// 	int n = 623;
// 	double a = k/2.;
// 	double b = 1. + (n-k)/2.;
// 	double res = my_gsl_cdf_beta_Pinv(alpha, a, b, 64, 1e-10);
//   printf("res: %g\n", res);
// 	res = my_gsl_cdf_beta_Pinv(alpha, a, b, 64, 1e-11);
// 	printf("res: %g\n", res);
// 	return 0;
// }
