// -----------------------------------------------------------------
// Print r_I (as defined in notes)
// corresponding to input integer displacement four-vector
// Use vegas to numerically calculate the infinite-volume Fourier transform
// Includes and definitions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"

#define    PI 3.141592653589793
#define TWOPI 6.283185307179586
#define  PISQ 9.869604401089358

// Double precision wall clock time in seconds
#include <sys/time.h>
double dclock() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}
// An alternative if sys/time.h isn't available
//double dclock() {
//  return (((double)clock()) / CLOCKS_PER_SEC);
//}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function to integrate -- only want real part of exp
// Delay sqrt until end -- negative evaluations seem possible
// Don't add any extra normalization (determinants cancel)
// Integrating over dp = dk / (2pi) removes 2pi factors from measure
static int f(const int *ndim, const cubareal *p,
             const int *ncomp, cubareal *ff, void *userdata) {

  int nx, ny, nz, nt;
  register cubareal tr, a, b, c, d;

  // Unpack userdata
  nx = ((int*)userdata)[0];
  ny = ((int*)userdata)[1];
  nz = ((int*)userdata)[2];
  nt = ((int*)userdata)[3];

  // Evaluate integrand
  tr = nx * p[0] + ny * p[1] + nz * p[2] + nt * p[3];
  *ff = cos(TWOPI * tr);

  a = sin(PI * p[0]);
  b = sin(PI * p[1]);
  c = sin(PI * p[2]);
  d = sin(PI * p[3]);
  tr = a * a + b * b + c * c + d * d;
  *ff /= tr;
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int i, n[4], nregions, fail;
  long long int neval;
  cubareal integral, error, prob, test, test_err;
  double abserr, dtime, val, err;

  // Read in command-line arguments
  // abserr is absolute error on the integral = 1 / r_I^2
  // It may need to be adjusted to obtain the desired r_I accuracy
  if (argc < 6) {
    printf("%s <nx> <ny> <nz> <nt> <error>\n", argv[0]);
    return -1;
  }
  for (i = 0; i < 4; i++)
    n[i] = atoi(argv[i + 1]);
  abserr = atof(argv[5]) / PISQ;

  // Strategy: Require consistency upon changing boundary from
  //           1e-8 to 1e-12.  If test passes, use latter result

  // Arguments: NDIM, NCOMP, function to integrate, user data, nvec
  //            EPSREL, EPSABS, VERBOSE+LEVEL, SEED
  //            MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
  //            BORDER, MAXCHISQ, MINDEVIATION,
  //            NGIVEN, LDXGIVEN=NDIM, XGIVEN, NEXTRA, PEAKFINDER
  //            STATEFILE, SPIN
  // KEY1=9 uses cubature rule of degree 9
  // KEY2=1 uses number of samples needed to reach desired accuracy
  // Points within BORDER of the integration boundary
  // will be evaluated by extrapolating from the interior
  dtime = -dclock();
  llDivonne(4, 1, f, (void*)n, 1,
            1e-8, abserr, 2, 0,
            0, 1e12, 9, 1, 1, 5,       // Need long long ints!
            1e-8, 10.0, 0.25,
            0, 4, NULL, 0, NULL,
            NULL, NULL,
            &nregions, &neval, &fail, &test, &test_err, &prob);
  llDivonne(4, 1, f, (void*)n, 1,
            1e-8, abserr, 2, 0,
            0, 1e12, 9, 1, 1, 5,       // Need long long ints!
            1e-12, 10.0, 0.25,
            0, 4, NULL, 0, NULL,
            NULL, NULL,
            &nregions, &neval, &fail, &integral, &error, &prob);

  val = fabs(test - integral);
  err = sqrt(test_err * test_err + error * error);
  if (val >= err) {
    printf("ERROR: Results sensitive to BORDER parameter:\n");
    printf("      |%.8g - %.8g| = %.8g\n", test, integral, val);
    printf("  vs. sqrt(%.4g^2 + %.4g^2) = %.4g\n", test_err, error, err);
    return -1;
  }
  printf("\nResults pass BORDER test:\n");
  printf("      |%.8g - %.8g| = %.8g\n", test, integral, val);
  printf("  vs. sqrt(%.4g^2 + %.4g^2) = %.4g\n", test_err, error, err);

  // Test has passed, so go ahead and print 1e-12 result
  integral *= PISQ;
  error *= PISQ;

  // Checked that 1 - prob matched confidence level for vegas
  if (fail == 0)
    printf("\nSuccess after %.2gM evalutions\n", (neval / 1e6));
  else
    printf("\nFailure after %2gM evalutions\n", (neval / 1e6));
  printf("result = %.8g %.4g with Q = %.2g\n",
         (double)integral, (double)error, 1 - (double)prob);

  // Print r_I itself with propagated uncertainty
  //  delta(1 / sqrt(r)) = delta(r) / 2r^(3 / 2)
  val = 1.0 / sqrt(integral);
  err = 0.5 * error * val * val * val;
  printf("--> r_I = %.8g %.4g\n", val, err);

  dtime += dclock();
  printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
