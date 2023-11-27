#include "./cloth_code.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double** table;

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

void initMatrix(int n, double UNUSED(mass), double UNUSED(fcon),
                int delta, double UNUSED(grav), double sep,
                double rball, double offset, double UNUSED(dt), double **x,
                double **y, double **z, double **cpx, double **cpy,
                double **cpz, double **fx, double **fy, double **fz,
                double **vx, double **vy, double **vz, double **oldfx,
                double **oldfy, double **oldfz) {
  int i, nx, ny;

  // Free any existing
  free(*x);
  free(*y);
  free(*z);
  free(*cpx);
  free(*cpy);
  free(*cpz);

  // allocate arrays to hold locations of nodes
  *x = (double *) malloc(n * n * sizeof(double));
  *y = (double *) malloc(n * n * sizeof(double));
  *z = (double *) malloc(n * n * sizeof(double));
  // This is for opengl stuff
  *cpx = (double *) malloc(n * n * sizeof(double));
  *cpy = (double *) malloc(n * n * sizeof(double));
  *cpz = (double *) malloc(n * n * sizeof(double));

  // initialize coordinates of cloth
  for (nx = 0; nx < n; nx++) {
    for (ny = 0; ny < n; ny++) {
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
  *fx = (double *) malloc(n * n * sizeof(double));
  *fy = (double *) malloc(n * n * sizeof(double));
  *fz = (double *) malloc(n * n * sizeof(double));
  *vx = (double *) malloc(n * n * sizeof(double));
  *vy = (double *) malloc(n * n * sizeof(double));
  *vz = (double *) malloc(n * n * sizeof(double));
  *oldfx = (double *) malloc(n * n * sizeof(double));
  *oldfy = (double *) malloc(n * n * sizeof(double));
  *oldfz = (double *) malloc(n * n * sizeof(double));
  for (i = 0; i < n * n; i++) {
    (*vx)[i] = 0.0;
    (*vy)[i] = 0.0;
    (*vz)[i] = 0.0;
    (*fx)[i] = 0.0;
    (*fy)[i] = 0.0;
    (*fz)[i] = 0.0;
  }
  // calculate table
  calculate_table(delta, sep);
}

void loopcode(int n, double mass, double fcon, int delta, double grav,
              double sep, double rball, double xball, double yball,
              double zball, double dt, double *x, double *y, double *z,
              double *fx, double *fy, double *fz, double *vx, double *vy,
              double *vz, double *oldfx, double *oldfy, double *oldfz,
              double *pe, double *ke, double *te) {
  int i, j, idx, n2 = n * n;
  double xdiff, ydiff, zdiff, vmag, damp, fvx, fvy, fvz, fmag, vfdot;
  double m_div = dt * 0.5 / mass, vmag_div, rball2 = rball * rball;

  // update position as per MD simulation
  for (idx = 0; idx < n2; idx++) {
      x[idx] += dt * (vx[idx] + fx[idx] * m_div);
      oldfx[idx] = fx[idx];
      y[idx] += dt * (vy[idx] + fy[idx] * m_div);
      oldfy[idx] = fy[idx];
      z[idx] += dt * (vz[idx] + fz[idx] * m_div);
      oldfz[idx] = fz[idx];

      xdiff = x[idx] - xball;
      ydiff = y[idx] - yball;
      zdiff = z[idx] - zball;
      
      vmag = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
      if (vmag < rball2) {
        vmag = sqrt(vmag);
        vmag_div = 1.0/vmag;
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


  *pe = eval_pef(n, delta, mass, grav, sep, fcon, x, y, z, fx, fy, fz);

  // Add a damping factor to eventually set velocity to zero
  damp = 0.995;
  *ke = 0.0;

  for (idx = 0; idx < n2; idx++) {
    vx[idx] = (vx[idx] +
                      (fx[idx] + oldfx[idx]) * m_div) *
                    damp;
    vy[idx] = (vy[idx] +
                      (fy[idx] + oldfy[idx]) * m_div) *
                    damp;
    vz[idx] = (vz[idx] +
                      (fz[idx] + oldfz[idx]) * m_div) *
                    damp;
    *ke += vx[idx] * vx[idx] + vy[idx] * vy[idx] +
            vz[idx] * vz[idx];

  }

  *ke = *ke * 0.5;
  *te = *pe + *ke;
}

double eval_pef(int n, int delta, double mass, double grav, double sep,
                double fcon, double *x, double *y, double *z, double *fx,
                double *fy, double *fz) 
{
  double pe, rlen, xdiff, ydiff, zdiff, vmag;
  double pe_delta, m_t_g = -mass * grav;
  int nx, ny, dx, dy, n2 = n * n, pef_idx, pef_idx2;

  pe = 0.0;
  // reset fx, fy, fz
  // use std::fill_n could be more efficient, e.g. std::fill_n(fx, n2, 0.0);
  for (pef_idx = 0; pef_idx < n2; pef_idx++) 
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
      for(ny = MAX(-dy, 0); ny < MIN(n - dy, n); ny++) 
      {
        for(nx = MAX(-dx, 0); nx < MIN(n - dx, n); nx++) 
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
