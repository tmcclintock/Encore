#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "allheader.h"

char SIMTYPE[1024];

FILE *fopen_retry(const char *filename, const char *mode)
{
  int try,Ntry = 10;
  FILE *fp;
  
  //try Ntry times, if opens, return fp
  for(try=0;try<Ntry;++try) {
    fp = fopen(filename,mode);
    
    if(fp != NULL)
      return fp;
  }
  
  //if we get to here, return NULL
  return NULL;
}

float read_period_length_file(char fname[])
{
  float (*read_pl)(char *);
  
  if(strcmp(SIMTYPE,"LGADGET") == 0)
    read_pl = &read_period_length_LGADGET;
  else {
    fprintf(stderr,"read period length - could not find I/O code for simulation type '%s'!\n",SIMTYPE);
    assert(0);
  }
  
  return read_pl(fname);
}

long read_num_chunks_file(char fname[])
{
  long (*read_nc)(char *);
  
  if(strcmp(SIMTYPE,"LGADGET") == 0)
    read_nc = &read_num_chunks_LGADGET;
  else {
    fprintf(stderr,"read num chunks - could not find I/O code for simulation type '%s'!\n",SIMTYPE);
    assert(0);
  }
  
  return read_nc(fname);
}

long read_num_parts_file(char fname[])
{
  long (*read_np)(char *);
  
  if(strcmp(SIMTYPE,"LGADGET") == 0)
    read_np = &read_num_parts_LGADGET;
  else {
    fprintf(stderr,"read num parts - could not find I/O code for simulation type '%s'!\n",SIMTYPE);
    assert(0);
  }
  
  return read_np(fname);
}

void read_parts_file(char fname[], float **px, float **py, float **pz, float **vx, float **vy , float **vz, long **id, long *Np)
{
  void (*read_file)(char *, float**, float**, float**, float**, float**, float**, long**, long*);
  
  if(strcmp(SIMTYPE,"LGADGET") == 0)
    read_file = &read_parts_LGADGET;
  else {
    fprintf(stderr,"read parts - could not find I/O code for simulation type '%s'!\n",SIMTYPE);
    assert(0);
  }
  
  read_file(fname,px,py,pz,vx,vy,vz,id,Np);
}

void write_parts_file(char fname_out[], char fname[], float *px, float *py, float *pz, float *vx, float *vy , float *vz, long *id, long Np)
{
  void (*write_file)(char*, char *, float*, float*, float*, float*, float*, float*, long*, long*);
  
  if(strcmp(SIMTYPE,"LGADGET") == 0)
    write_file = &write_parts_LGADGET;
  else {
    fprintf(stderr,"write parts - could not find I/O code for simulation type '%s'!\n",SIMTYPE);
    assert(0);
  }
  
  write_file(fname_out,fname,px,py,pz,vx,vy,vz,id,Np);
}

