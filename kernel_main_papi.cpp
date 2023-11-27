#ifdef PAPI
  #include <papi.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./cloth_code.h"
#include "./cloth_param.h"

#define NUM 9
long long counters[NUM], start, stop;
int EventSet = PAPI_NULL;
int Events[NUM] = {
    PAPI_TOT_INS, //0
    PAPI_TOT_CYC, //1
    PAPI_L1_DCM,  //2
    PAPI_L3_TCA,  //3
    PAPI_L3_TCM,  //4
    PAPI_L2_DCM,  //5
    PAPI_L2_DCA,  //6
    PAPI_BR_MSP,  //7
    PAPI_BR_CN,   //8
};

int main(int argc, char **argv) {
  int i, iter;
  double pe, ke, te;
  // assess input flags

  for (i = 1; i < argc; i += 2) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'n':
        n = atoi(argv[i + 1]);
        break;
      case 's':
        sep = atof(argv[i + 1]);
        break;
      case 'm':
        mass = atof(argv[i + 1]);
        break;
      case 'f':
        fcon = atof(argv[i + 1]);
        break;
      case 'd':
        delta = atoi(argv[i + 1]);
        break;
      case 'g':
        grav = atof(argv[i + 1]);
        break;
      case 'b':
        rball = atof(argv[i + 1]);
        break;
      case 'o':
        offset = atof(argv[i + 1]);
        break;
      case 't':
        dt = atof(argv[i + 1]);
        break;
      case 'i':
        maxiter = atoi(argv[i + 1]);
        break;
      default:
        printf(" %s\n"
               "Nodes_per_dimension:             -n int \n"
               "Grid_separation:                 -s float \n"
               "Mass_of_node:                    -m float \n"
               "Force_constant:                  -f float \n"
               "Node_interaction_level:          -d int \n"
               "Gravity:                         -g float \n"
               "Radius_of_ball:                  -b float \n"
               "offset_of_falling_cloth:         -o float \n"
               "timestep:                        -t float \n"
               "num iterations:                  -i int \n",
               argv[0]);
        return -1;
      }
    } else {
      printf(" %s\n"
             "Nodes_per_dimension:             -n int \n"
             "Grid_separation:                 -s float \n"
             "Mass_of_node:                    -m float \n"
             "Force_constant:                  -f float \n"
             "Node_interaction_level:          -d int \n"
             "Gravity:                         -g float \n"
             "Radius_of_ball:                  -b float \n"
             "offset_of_falling_cloth:         -o float \n"
             "timestep:                        -t float \n"
             "num iterations:                  -i int \n",
             argv[0]);
      return -1;
    }
  }

  initMatrix(n, mass, fcon, delta, grav, sep, rball, offset, dt, &x, &y, &z,
             &cpx, &cpy, &cpz, &fx, &fy, &fz, &vx, &vy, &vz, &oldfx, &oldfy,
             &oldfz);

  // papi
  start = PAPI_get_virt_usec();
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_create_eventset(&EventSet);
  PAPI_add_events(EventSet, Events, NUM);
  PAPI_start(EventSet);

  for (iter = 0; iter < maxiter; iter++) {
    loopcode(n, mass, fcon, delta, grav, sep, rball, xball, yball, zball, dt, x,
             y, z, fx, fy, fz, vx, vy, vz, oldfx, oldfy, oldfz, &pe, &ke, &te);
  }
  stop = PAPI_get_virt_usec();
  PAPI_stop(EventSet, counters);
  printf("%d,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%d,%lld,", n, sep, mass, fcon, delta, grav, rball, offset, dt, maxiter, stop - start);
  for (int i = 0; i < NUM; i++) {
    printf("%lld", counters[i]);
    if (i != NUM - 1) printf(",");
  }
  printf("\n");

  return 0;
}
