#include "./cloth_code_omp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __SSE2__
#include <immintrin.h>
#define _mul _mm256_mul_pd
#define _div _mm256_div_pd
#define _add _mm256_add_pd
#define _sub _mm256_sub_pd
#define _set1 _mm256_set1_pd
#define _load _mm256_loadu_pd
#define _store _mm256_storeu_pd
#endif
double** table;

double sum_elements(__m256d vec)
{
  // Horizontal add
  __m256d sum = _mm256_hadd_pd(vec, vec); // This produces [a0+a1, a2+a3, a0+a1, a2+a3]

  // Extract upper 128 bits of the result
  __m128d high = _mm256_extractf128_pd(sum, 1);
  __m128d low = _mm256_castpd256_pd128(sum);
  __m128d finalSum = _mm_add_pd(high, low);

  double result[2];
  _mm_storeu_pd(result, finalSum);

  return result[0]; 
}

void calculate_table(int delta, double sep) {
  int ii, jj;
  // malloc for table
  table = (double**)malloc((delta + 1) * sizeof(double *));
  for(ii = 0; ii < delta + 1; ii++) {
    table[ii] = (double*)malloc((delta + 1) * sizeof(double));
  }
  for (ii = 0; ii < delta + 1; ii++) {
    for (jj = ii; jj < delta + 1; jj++) {
      table[ii][jj] = table[jj][ii] = sqrt((double)(ii * ii + jj * jj)) * sep;
    }
  }
}

double query_table(int a, int b) {
  a = ((a >> 31) ^ a) - (a >> 31);
  b = ((b >> 31) ^ b) - (b >> 31);
  return table[a][b];
}

__m256d ONE = _set1(1.0), ZERO = _mm256_setzero_pd(), ZPO = _set1(0.1);
void initMatrix(int n, double UNUSED(mass), double UNUSED(fcon),
                int delta, double UNUSED(grav), double sep,
                double rball, double offset, double UNUSED(dt), double **x,
                double **y, double **z, double **cpx, double **cpy,
                double **cpz, double **fx, double **fy, double **fz,
                double **vx, double **vy, double **vz, double **oldfx,
                double **oldfy, double **oldfz, int tn)
{
  int i, nx, ny;

  // Free any existing
  free(*x);
  free(*y);
  free(*z);
  free(*cpx);
  free(*cpy);
  free(*cpz);

  // allocate arrays to hold locations of nodes
  *x = (double *)malloc(n * n * sizeof(double));
  *y = (double *)malloc(n * n * sizeof(double));
  *z = (double *)malloc(n * n * sizeof(double));
  // This is for opengl stuff
  *cpx = (double *)malloc(n * n * sizeof(double));
  *cpy = (double *)malloc(n * n * sizeof(double));
  *cpz = (double *)malloc(n * n * sizeof(double));

  // initialize coordinates of cloth
  for (nx = 0; nx < n; nx++)
  {
    for (ny = 0; ny < n; ny++)
    {
      (*x)[n * nx + ny] = nx * sep - (n - 1) * sep * 0.5 + offset;
      (*z)[n * nx + ny] = rball + 1;
      (*y)[n * nx + ny] = ny * sep - (n - 1) * sep * 0.5 + offset;
      (*cpx)[n * nx + ny] = 0;
      (*cpz)[n * nx + ny] = 1;
      (*cpy)[n * nx + ny] = 0;
    }
  }

  // Throw away existing arrays
  free(*fx);
  free(*fy);
  free(*fz);
  free(*vx);
  free(*vy);
  free(*vz);
  free(*oldfx);
  free(*oldfy);
  free(*oldfz);
  // Alloc new
  *fx = (double *)malloc(n * n * sizeof(double));
  *fy = (double *)malloc(n * n * sizeof(double));
  *fz = (double *)malloc(n * n * sizeof(double));
  *vx = (double *)malloc(n * n * sizeof(double));
  *vy = (double *)malloc(n * n * sizeof(double));
  *vz = (double *)malloc(n * n * sizeof(double));
  *oldfx = (double *)malloc(n * n * sizeof(double));
  *oldfy = (double *)malloc(n * n * sizeof(double));
  *oldfz = (double *)malloc(n * n * sizeof(double));
  for (i = 0; i < n * n; i++)
  {
    (*vx)[i] = 0.0;
    (*vy)[i] = 0.0;
    (*vz)[i] = 0.0;
    (*fx)[i] = 0.0;
    (*fy)[i] = 0.0;
    (*fz)[i] = 0.0;
  }
  calculate_table(delta, sep);
}

void loopcode(int n, double mass, double fcon, int delta, double grav,
              double sep, double rball, double xball, double yball,
              double zball, double dt, double *x, double *y, double *z,
              double *fx, double *fy, double *fz, double *vx, double *vy,
              double *vz, double *oldfx, double *oldfy, double *oldfz,
              double *pe, double *ke, double *te, int tn)
{
  int i, j, idx, n2 = n * n;
  double xdiff, ydiff, zdiff, vmag, damp, fvx, fvy, fvz, fmag, vfdot;
  double m_div = dt * 0.5 / mass, vmag_div, rball2 = rball * rball;
  __m256d _fx, _fy, _fz, _vx, _vy, _vz, _x, _y, _z, _vmag, _xdiff, _ydiff, _zdiff, _fvx, _fvy, _fvz, _fmag, _vfdot, _ke_delta;
  __m256d _rball = _set1(rball), _rball2 = _set1(rball * rball), _dt = _set1(dt), _xball = _set1(xball), _yball = _set1(yball), _zball = _set1(zball), _m_div = _set1(m_div), _damp = _set1(0.995), _inv_2mass = _set1(m_div);
  __m256d _cmp_res, _cmp_res1;

  #pragma omp parallel for lastprivate(idx) private(_fx, _fy, _fz, _vx, _vy, _vz, _x, _y, _z, _xdiff, _ydiff, _zdiff, _vmag, _cmp_res, _fvx, _fvy, _fvz, _fmag, _vfdot) num_threads(tn) schedule(static)
  for (idx = 0; idx <= n2 - 4; idx += 4)
  {
    _fx = _load(fx + idx);
    _fy = _load(fy + idx);
    _fz = _load(fz + idx);

    _vx = _load(vx + idx);
    _vy = _load(vy + idx);
    _vz = _load(vz + idx);

    _x = _add(_load(x + idx), _mul(_dt, _add(_mul(_fx, _m_div), _vx)));
    _y = _add(_load(y + idx), _mul(_dt, _add(_mul(_fy, _m_div), _vy)));
    _z = _add(_load(z + idx), _mul(_dt, _add(_mul(_fz, _m_div), _vz)));
    _store(x + idx, _x);
    _store(y + idx, _y);
    _store(z + idx, _z);

    _store(oldfx + idx, _fx);
    _store(oldfy + idx, _fy);
    _store(oldfz + idx, _fz);

    _xdiff = _sub(_x, _xball);
    _ydiff = _sub(_y, _yball);
    _zdiff = _sub(_z, _zball);
    _vmag = _add(_add(_mul(_xdiff, _xdiff), _mul(_ydiff, _ydiff)), _mul(_zdiff, _zdiff));
    _cmp_res = _mm256_cmp_pd(_vmag, _rball2, _CMP_LT_OQ);
    
    if (_mm256_movemask_pd(_cmp_res) == 0)
    {
      continue;
    }
    _vmag = _mm256_sqrt_pd(_vmag);
    _vmag = _mm256_div_pd(ONE, _vmag);

    _fvx = _mul(_mul(_xdiff, _rball), _vmag);
    _fvy = _mul(_mul(_ydiff, _rball), _vmag);
    _fvz = _mul(_mul(_zdiff, _rball), _vmag);

    _fmag = _add(_add(_mul(_fvx, _fvx), _mul(_fvy, _fvy)), _mul(_fvz, _fvz));

    _store(x + idx, _mm256_blendv_pd(_x, _add(_xball, _fvx), _cmp_res));
    _store(y + idx, _mm256_blendv_pd(_y, _add(_yball, _fvy), _cmp_res));
    _store(z + idx, _mm256_blendv_pd(_z, _add(_zball, _fvz), _cmp_res));
    _vfdot = _div(_add(_add(_mul(_fvx, _vx), _mul(_fvy, _vy)), _mul(_fvz, _vz)), _fmag);
    _store(vx + idx, _mm256_blendv_pd(_vx, _mul(ZPO, _sub(_vx, _mul(_vfdot, _fvx))), _cmp_res));
    _store(vy + idx, _mm256_blendv_pd(_vy, _mul(ZPO, _sub(_vy, _mul(_vfdot, _fvy))), _cmp_res));
    _store(vz + idx, _mm256_blendv_pd(_vz, _mul(ZPO, _sub(_vz, _mul(_vfdot, _fvz))), _cmp_res));
  }
  // #pragma omp parallel for private(idx, xdiff, ydiff, zdiff, vmag, vmag_div, fvx, fvy, fvz, fmag, vfdot) num_threads(tn) schedule(static)
  for (; idx < n2; idx++)
  {
    // update position as per MD simulation
    x[idx] += dt * (vx[idx] + fx[idx] * m_div);
    oldfx[idx] = fx[idx];
    y[idx] += dt * (vy[idx] + fy[idx] * m_div);
    oldfy[idx] = fy[idx];
    z[idx] += dt * (vz[idx] + fz[idx] * m_div);
    oldfz[idx] = fz[idx];
    //	apply constraints - push cloth outside of ball
    xdiff = x[idx] - xball;
    ydiff = y[idx] - yball;
    zdiff = z[idx] - zball;

    vmag = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
    if (vmag < rball2)
    {
      vmag = sqrt(vmag);
      vmag_div = 1.0 / vmag;
      fvx = xdiff * rball * vmag_div;
      fvy = ydiff * rball * vmag_div;
      fvz = zdiff * rball * vmag_div;
      // fmag = sqrt(fvx * fvx + fvy * fvy + fvz * fvz);
      fmag = fvx * fvx + fvy * fvy + fvz * fvz;
      x[idx] = xball + fvx;
      y[idx] = yball + fvy;
      z[idx] = zball + fvz;
      // update velocity
      // vfdot = (fvx * vx[idx] + fvy * vy[idx] + fvz * vz[idx]) / fmag / fmag;
      vfdot = (fvx * vx[idx] + fvy * vy[idx] + fvz * vz[idx]) / fmag;
      vx[idx] = 0.1 * (vx[idx] - vfdot * fvx);
      vy[idx] = 0.1 * (vy[idx] - vfdot * fvy);
      vz[idx] = 0.1 * (vz[idx] - vfdot * fvz);
    }
  }
  // }

  *pe = eval_pef(n, delta, mass, grav, sep, fcon, x, y, z, fx, fy, fz, tn);

  // Add a damping factor to eventually set velocity to zero
  damp = 0.995;
  *ke = 0.0;
  double kke = 0.0;
  #pragma omp parallel for lastprivate(idx) reduction(+:kke) private(_vx, _vy, _vz, _ke_delta) num_threads(tn) schedule(static)
  for (idx = 0; idx <= n2 - 4; idx += 4)
  {
    _store(vx + idx, _mul(_damp, _add(_load(vx + idx), _mul(_add(_load(fx + idx), _load(oldfx + idx)), _inv_2mass))));
    _store(vy + idx, _mul(_damp, _add(_load(vy + idx), _mul(_add(_load(fy + idx), _load(oldfy + idx)), _inv_2mass))));
    _store(vz + idx, _mul(_damp, _add(_load(vz + idx), _mul(_add(_load(fz + idx), _load(oldfz + idx)), _inv_2mass))));

    _vx = _load(vx + idx);
    _vy = _load(vy + idx);
    _vz = _load(vz + idx);
    _ke_delta = _add(_add(_mul(_vx, _vx), _mul(_vy, _vy)), _mul(_vz, _vz));
    kke += sum_elements(_ke_delta);
  }

  // #pragma omp parallel for private(idx) reduction(+:kke)
  // #pragma omp simd reduction(+:kke) aligned(vx: 64, vy: 64, vz: 64, fx: 64, fy: 64, fz: 64, oldfx: 64, oldfy: 64, oldfz: 64)
  for (; idx < n2; idx++)
  {
    vx[idx] = (vx[idx] +
               (fx[idx] + oldfx[idx]) * m_div) *
              damp;
    vy[idx] = (vy[idx] +
               (fy[idx] + oldfy[idx]) * m_div) *
              damp;
    vz[idx] = (vz[idx] +
               (fz[idx] + oldfz[idx]) * m_div) *
              damp;

    kke += vx[idx] * vx[idx] + vy[idx] * vy[idx] +
           vz[idx] * vz[idx];
  }
  // }
  *ke = kke * 0.5;
  *te = *pe + *ke;
}

double eval_pef(int n, int delta, double mass, double grav, double sep,
                double fcon, double *x, double *y, double *z, double *fx,
                double *fy, double *fz, int tn)
{
  double pe, rlen, xdiff, ydiff, zdiff, vmag;
  double pe_delta, m_t_g = -mass * grav;

  __m256d _xdiff, _ydiff, _zdiff, _vmag, _fx, _fy, _fz, _nx, _ny, _x, _y, _z, _rlen, _dx, _dy, _pe_delta, _pe_delta1, _v_r;
  __m256d _m_t_g = _set1(m_t_g), _sep = _set1(sep), _fcon = _set1(fcon);
  __m256d _cmp_res;

  int nx, ny, dx, dy, start_ny, start_nx, end_nx, end_ny, pef_idx, pef_idx2, n2 = n * n;
  pe = 0.0;
  
  // reset fx, fy, fz
  // use std::fill_n could be more efficient, e.g. std::fill_n(fx, n2, 0.0);
  for (pef_idx = 0; pef_idx < n2 - 4; pef_idx+=4) 
  {
    _store(fx + pef_idx, ZERO);
    _store(fy + pef_idx, ZERO);
    _store(fz + pef_idx, _m_t_g);
  }
  #pragma omp simd aligned(fx:64, fy: 64, fz: 64)
  for (; pef_idx < n2; pef_idx++) 
  {
    fx[pef_idx] = 0.0;
    fy[pef_idx] = 0.0;
    fz[pef_idx] = m_t_g;
  }
  
  for (dy = -delta; dy <= delta; dy++) 
  {
    for (dx = -delta; dx <= delta; dx++) 
    {
      if (dy == 0 && dx == 0) {
        continue;
      }
      rlen = query_table(dx, dy);
      start_nx = MAX(-dx, 0);
      start_ny = MAX(-dy, 0);
      end_nx = MIN(n - dx, n);
      end_ny = MIN(n - dy, n);
      _rlen = _set1(rlen);
      #pragma omp parallel for lastprivate(nx) private(pef_idx, pef_idx2, _xdiff, _ydiff, _zdiff, _vmag, _v_r, _pe_delta, _pe_delta1) reduction(+:pe) num_threads(tn) collapse(2) schedule(static)
      for(ny = start_ny; ny < end_ny; ny++) 
      {
        for(nx = start_nx; nx <= end_nx - 4; nx+=4) 
        {
          pef_idx = ny * n + nx;
          pef_idx2 = (dy + ny) * n + (dx + nx);
          // compute actual distance
          _xdiff = _sub(_load(x + pef_idx2), _load(x + pef_idx));
          _ydiff = _sub(_load(y + pef_idx2), _load(y + pef_idx));
          _zdiff = _sub(_load(z + pef_idx2), _load(z + pef_idx));
          _vmag = _mm256_sqrt_pd(_add(_add(_mul(_xdiff, _xdiff), _mul(_ydiff, _ydiff)), _mul(_zdiff, _zdiff)));

          _v_r = _sub(_vmag, _rlen);
          _pe_delta = _mul(_fcon, _mul(_v_r, _v_r));
          _pe_delta1 = _mul(_fcon, _div(_v_r, _vmag));

          pe += sum_elements(_pe_delta);
          // _xdiff = ;
          // _ydiff = ;
          // _zdiff = ;
          _store(fx + pef_idx, _add(_load(fx + pef_idx), _mul(_xdiff, _pe_delta1)));
          _store(fy + pef_idx, _add(_load(fy + pef_idx), _mul(_ydiff, _pe_delta1)));
          _store(fz + pef_idx, _add(_load(fz + pef_idx), _mul(_zdiff, _pe_delta1)));
        }
        
      }
      #pragma omp simd private(pef_idx, pef_idx2, nx, ny, xdiff, ydiff, zdiff, vmag, pe_delta) reduction(+:pe, fx[:n2], fy[:n2], fz[:n2]) collapse(2)
      for(ny = start_ny; ny < end_ny; ny++) 
      {
        for(nx = end_nx - ((end_nx - start_nx) % 4); nx < end_nx; nx++) 
        {
          pef_idx = ny * n + nx;
          pef_idx2 = (dy + ny) * n + (dx + nx);
          // compute actual distance
          xdiff = x[pef_idx2] - x[pef_idx];
          ydiff = y[pef_idx2] - y[pef_idx];
          zdiff = z[pef_idx2] - z[pef_idx];
          vmag = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
          // potential energy and force
          pe_delta = fcon * (vmag - rlen);
          pe += pe_delta * (vmag - rlen);
          // update fx
          pe_delta /= vmag;
          fx[pef_idx] += xdiff * pe_delta;
          fy[pef_idx] += ydiff * pe_delta;
          fz[pef_idx] += zdiff * pe_delta;
        }
      }
    }
  }
  return 0.5 * pe;
}
