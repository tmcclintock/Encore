//1400^3 = 2 744 000 000 particles total!!

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "allheader.h"

int main(int argc, char **argv)
{
  if(argc != 5) {
    fprintf(stderr,"subsamp_parts - code to sub-sample snapshots\nusage: subsamp_parts SIMTYPE /path/to/files/filebase NumToKeep /path/to/output\n\
Supported simulation types are\n\
    LGADGET\n");
    exit(1);
  }
  
  assert(argc == 5);
  sprintf(SIMTYPE,"%s",argv[1]);
  char *fbase = argv[2];
  long num_keep = atol(argv[3]);
  char fname[MAX_FILENAME];
  char *fname_out = argv[4];
  
  double frac;
  long num_parts;
  long num_have = 0;
  float *px,*py,*pz,*vx,*vy,*vz;
  long *id;
  
  long chunk,num_chunks,i;
  float *pxi,*pyi,*pzi,*vxi,*vyi,*vzi;
  long *idi;
  long Npi;
  
  //do RNG
  const gsl_rng_type *T;
  gsl_rng *r;  
  T = gsl_rng_ranlxd1;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r,(unsigned long int) (time(NULL)));
  
  //get file props and fracs
  sprintf(fname,"%s.0",fbase);
  num_chunks = read_num_chunks_file(fname);
  num_parts = read_num_parts_file(fname);
  
  if(num_parts < num_keep) {
    num_keep = num_parts;
    frac = 1.1;
  } else {
    frac = ((double) num_keep)/((double) num_parts);
  }
  
  //io info
  fprintf(stderr,"fbase %s in %ld chunks\n",fbase,num_chunks);
  fprintf(stderr,"sub samp %g, num_parts = %ld, num_keep = %ld\n",frac,num_parts,num_keep);
  fflush(stderr);
  
  //alloc mem
  px = (float*)malloc(sizeof(float)*num_keep);
  assert(px != NULL);
  py = (float*)malloc(sizeof(float)*num_keep);
  assert(py != NULL);
  pz = (float*)malloc(sizeof(float)*num_keep);
  assert(pz != NULL);  
  vx = (float*)malloc(sizeof(float)*num_keep);
  assert(vx != NULL);
  vy = (float*)malloc(sizeof(float)*num_keep);
  assert(vy != NULL);
  vz = (float*)malloc(sizeof(float)*num_keep);
  assert(vz != NULL);  
  id = (long*)malloc(sizeof(long)*num_keep);
  
  //read each chunk, then keep right particles
  for(chunk=0;chunk<num_chunks;++chunk) {
    fprintf(stderr,"chunk %ld of %ld: %ld parts\n",chunk+1,num_chunks,num_have);
    fflush(stderr);
    
    sprintf(fname,"%s.%ld",fbase,chunk);    
    read_parts_file(fname,&pxi,&pyi,&pzi,&vxi,&vyi,&vzi,&idi,&Npi);
    
    for(i=0;i<Npi;++i) {
      if(gsl_rng_uniform(r) < frac) {
	px[num_have] = pxi[i];
	py[num_have] = pyi[i];
	pz[num_have] = pzi[i];
	
	vx[num_have] = vxi[i];
	vy[num_have] = vyi[i];
	vz[num_have] = vzi[i];

	id[num_have] = idi[i];
	
	++num_have;
	
	if(num_have == num_keep)
	  break;
      }
    }
    
    free(pxi);
    free(pyi);
    free(pzi);
    free(vxi);
    free(vyi);
    free(vzi);
    free(idi);
    
    if(num_have == num_keep)
      break;
  }
  
  //write data
  sprintf(fname,"%s.0",fbase);    
  write_parts_file(fname_out,fname,px,py,pz,vx,vy,vz,id,num_have);
  
  //clean up
  gsl_rng_free(r);
  free(px);
  free(py);
  free(pz);
  free(vx);
  free(vy);
  free(vz);
  free(id);
  
  return 0;
}

