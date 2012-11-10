#include <math.h>
#include <mex.h>

// supply your own random number!
// k = fast_sample(pdf, r)
// pdf and r must be double, and one-dimensional.
// No error checking!
#define mxGetSize(m)    (mxGetN(m) * mxGetM(m))
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  if (nlhs != 1 || nrhs != 2) {
    mexErrMsgTxt("Usage: k = fast_sample(pdf, r). No error checking!");
  }
  size_t sz   = mxGetSize(prhs[0]);
  double *pdf = mxGetPr(prhs[0]);
  double r    = mxGetScalar(prhs[1]);

  register double sum = 0.0;
  double cdf[sz];
  for (int i = 0; i < sz; i++) {
    sum += pdf[i];
    cdf[i] = sum;
  }

  double target = sum * r;
  for (int i = 0; i < sz; i++) {
    if (cdf[i] > target) {
      // Remember to + 1 for 1-based index.
      plhs[0] = mxCreateDoubleScalar(i + 1);
      return;
    }
  }
  plhs[0] = mxCreateDoubleScalar(sz);
}
