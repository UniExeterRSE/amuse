#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include "query_struct.h"

typedef struct payload_t
{
    struct sph_particle_data sph;
    struct particle_data p;
} payload_t;


typedef struct data_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

// Input arrays
static double *point_req_x = NULL;
static double *point_req_y = NULL;
static double *point_req_z = NULL;
static double *nearest_hsml = NULL;


typedef struct data_out
{
  payload_t payload;
  double Distance;
} data_out;

// Output arrays
static payload_t *payload_responses = NULL;
static double *nearest_dist = NULL;


static void particle2in(data_in * in, int place, int firstnode)
{
  in->Pos[0] = point_req_x[place];
  in->Pos[1] = point_req_y[place];
  in->Pos[2] = point_req_z[place];
  in->Hsml = nearest_hsml[place];
  in->Firstnode = firstnode;
}

static void unpack(const data_out * const out, int place)
{
  nearest_dist[place] = out->Distance;
  payload_responses[place] = out->payload;
}


static void out2particle(const data_out * const out, int place, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)  // initial store
    {
      unpack(out, place);
    }
  else  // combine
    {
      if(out->Distance < nearest_dist[place])
        unpack(out, place);
    }
}



static data_in *DataIn = NULL;
static data_in *DataGet = NULL;
static data_out *DataResult = NULL;
static data_out *DataOut = NULL;

#include "src/utils/generic_comm_helpers2.h"

static int Ncount = 0;

static int find_nearest_evaluate(int target, int mode, int thread_id);

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Ncount))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ncount)
          break;

        if(nearest_dist[i] > 1.0e29)
          find_nearest_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        find_nearest_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


static double obtain_fields_at_points(int npoints);


void query_struct(const int npoints, const double x[], const double y[], const double z[],
    struct sph_particle_data sph_out[], struct particle_data p_out[])
{
  // Allocate request and response arrays
  point_req_x = (double *) mymalloc("query_struct_req_x", sizeof(double) * npoints);
  point_req_y = (double *) mymalloc("query_struct_req_y", sizeof(double) * npoints);
  point_req_z = (double *) mymalloc("query_struct_req_z", sizeof(double) * npoints);

  payload_responses = (payload_t *) mymalloc("payload_response", sizeof(payload_t) * npoints);

  // Copy requests to input arrays
  for (int i=0; i<npoints; i++) {
    point_req_x[i] = x[i];
    point_req_y[i] = y[i];
    point_req_z[i] = z[i];
  }

  // Obtain fields at points
  const double t_exch = obtain_fields_at_points(npoints);

  // Copy back data to output arrays
  for (int i=0; i<npoints; i++) {
    sph_out[i] = payload_responses[i].sph;
    p_out[i] = payload_responses[i].p;
  }

  myfree(payload_responses); payload_responses = NULL;

  myfree(point_req_z); point_req_z = NULL;
  myfree(point_req_y); point_req_y = NULL;
  myfree(point_req_x); point_req_x = NULL;
}


void query_density(const int npoints, const double x[], const double y[], const double z[], double density[])
{
    auto sph = new struct sph_particle_data[npoints];
    auto p = new struct particle_data[npoints];
    query_struct(npoints, x, y, z, sph, p);
    for (int i=0; i<npoints; i++) {
        density[i] = sph[i].Density;
    }
    delete [] sph;
    delete [] p;
}

static double obtain_fields_at_points(int npoints)
{
  double tstart = second();
  long long ntot, npleft, i;
  int iter;

  mpi_printf("start collecting data (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  nearest_dist = (double *) mymalloc("query_nearest_dist", sizeof(double) * npoints);
  nearest_hsml = (double *) mymalloc("query_nearest_hsml", sizeof(double) * npoints);

  for(int n=0; n<npoints; n++)
    {
      nearest_dist[n] = 1.0e30;
      nearest_hsml[n] = All.BoxSize / pow(All.TotNumGas, 1.0 / 3);
    }

  Ncount = npoints;

  generic_set_MaxNexport();

  iter = 0;
  do {
  generic_comm_pattern(Ncount, kernel_local, kernel_imported);

  /* do final operations on results */
  for(i = 0, npleft = 0; i < Ncount; i++)
    {
      if(nearest_dist[i] > 1.0e29)
        {
          /* need to redo this particle */
          npleft++;
          nearest_hsml[i] *= 2.0;
        }
       else
        {
          nearest_dist[i] = 0;    /* we not continue to search for this particle */
        }
     }
    sumup_longs(1, &npleft, &ntot);
    if(ntot > 0)
      {
        iter++;
        if(iter > 0 && ThisTask == 0)
          {
            printf("query_struct: nearest iteration %d: need to repeat for %lld particles.\n", iter, ntot);
            myflush(stdout);
          }

        if(iter > MAXITER)
          terminate("query_struct: failed to converge");
      }
    }
  while(ntot > 0);

  myfree(nearest_hsml);
  myfree(nearest_dist);

  double tend = second();
  return timediff(tstart, tend);
}


static int find_nearest_evaluate(int target, int mode, int thread_id)
{
  int j, n;
  int numnodes, *firstnode;

  MyFloat h;
  double r2max;
  double dx, dy, dz, r2;
  MyDouble *pos;

  data_in local, *target_data;
  data_out out;


  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  h = target_data->Hsml;

  int index = -1;
  r2max = 1.0e60;
  out.Distance = sqrt(r2max);

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  //int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int mode, int thread_id, int numnodes,
                                  //int *firstnode);

  if(nfound < 0)
    return -1;

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];
      double xtmp, ytmp, ztmp;
      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < r2max && r2 < h * h)
        {
          index = j;
          r2max = r2;
        }
    }

  if(index >= 0)
    {
      if(index >= NumGas)
        terminate("index >= NumGas");

      out.Distance = sqrt(r2max);

      out.payload.sph = SphP[index];
      out.payload.p = P[index];
    }
  else
    {
      out.Distance = 1.e60;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}
