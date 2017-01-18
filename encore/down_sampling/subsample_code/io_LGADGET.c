#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "allheader.h"

//LGadget-2 I/O header
struct io_header_1 {
  unsigned int npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int flag_feedback;  /*!< flags whether feedback from star formation is included */
  unsigned int npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores                                                 
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int flag_cooling;   /*!< flags whether radiative cooling is included */
  int num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  int flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int flag_metals;         /*!< flags whether metal enrichment is included */
  int hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  //char fill[84];              /*!< fills to 256 Bytes */                                                                                                                         
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of par                                                                                                          
                                         ticles of each type */
  char fill[60];
};

float read_period_length_LGADGET(char fname[])
{
  FILE *fp;
  struct io_header_1 header1;
  int dummy;
#define SKIP fread(&dummy,sizeof(dummy),(size_t) 1,fp);

  fp = fopen_retry(fname,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname);
    assert(fp != NULL);
  }
  
  SKIP;
  fread(&header1,sizeof(header1),(size_t) 1,fp);
  SKIP;
  
  fclose(fp);

#undef SKIP

  return header1.BoxSize;
}

long read_num_parts_LGADGET(char fname[])
{
  FILE *fp;
  struct io_header_1 header1;
  int dummy;
#define SKIP fread(&dummy,sizeof(dummy),(size_t) 1,fp);

  fp = fopen_retry(fname,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname);
    assert(fp != NULL);
  }
  
  SKIP;
  fread(&header1,sizeof(header1),(size_t) 1,fp);
  SKIP;
  
  fclose(fp);
  
#undef SKIP
  
  return (long) (header1.npartTotal[1]) + ((long) (header1.npartTotal[2]) << 32);
}

long read_num_chunks_LGADGET(char fname[])
{
  FILE *fp;
  struct io_header_1 header1;
  int dummy;
#define SKIP fread(&dummy,sizeof(dummy),(size_t) 1,fp);

  fp = fopen_retry(fname,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname);
    assert(fp != NULL);
  }
  
  SKIP;
  fread(&header1,sizeof(header1),(size_t) 1,fp);
  SKIP;
  printf("Redshift = %f\n",header1.redshift);
  
  fclose(fp);
  
#undef SKIP
  
  return (long) (header1.num_files);
}

void read_parts_LGADGET(char fname[], float **px, float **py, float **pz, float **vx, float **vy , float **vz, long **id, long *Np)
{
  //vars and macros
  FILE *fp;
  struct io_header_1 header1;
  int dummy_start,dummy_stop,k;
  float *partChunk;
  long *idChunk_long;
  int *idChunk_int;

#define START fread(&dummy_start,sizeof(dummy_start),(size_t) 1,fp);
#define STOP fread(&dummy_stop,sizeof(dummy_stop),(size_t) 1,fp);
#define CHECK(SZE,NP) assert(dummy_start == dummy_stop && (SZE)*(NP) == dummy_start);
  
  //open file
  fp = fopen_retry(fname,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname);
    assert(fp != NULL);
  }
  
  //do header
  START;
  fread(&header1,sizeof(header1),(size_t) 1,fp);
  STOP;
  CHECK(256,1);
  
  //alloc mem
  *Np = header1.npart[1];

  *px = (float*)malloc(sizeof(float)*(*Np));
  assert((*px) != NULL);
  *py = (float*)malloc(sizeof(float)*(*Np));
  assert((*py) != NULL);
  *pz = (float*)malloc(sizeof(float)*(*Np));
  assert((*pz) != NULL);
  
  *vx = (float*)malloc(sizeof(float)*(*Np));
  assert((*vx) != NULL);
  *vy = (float*)malloc(sizeof(float)*(*Np));
  assert((*vy) != NULL);
  *vz = (float*)malloc(sizeof(float)*(*Np));
  assert((*vz) != NULL);
  
  *id = (long*)malloc(sizeof(long)*(*Np));

  partChunk = (float*)malloc(sizeof(float)*(*Np)*3);
  assert(partChunk != NULL);
  idChunk_long = (long*)partChunk;
  idChunk_int = (int*)partChunk;
  
  //read positions                                                                                                                                                                                           
  START;
  assert(fread(partChunk,(size_t) (*Np),3*sizeof(float),fp) == 3*sizeof(float));
  for(k=0;k<(*Np);++k) {
    (*px)[k] = partChunk[k*3 + 0];
    (*py)[k] = partChunk[k*3 + 1];
    (*pz)[k] = partChunk[k*3 + 2];
    
    while((*px)[k] < 0.0)
      (*px)[k] += header1.BoxSize;
    while((*px)[k] > header1.BoxSize)
      (*px)[k] -= header1.BoxSize;
    
    while((*py)[k] < 0.0)
      (*py)[k] += header1.BoxSize;
    while((*py)[k] > header1.BoxSize)
      (*py)[k] -= header1.BoxSize;
    
    while((*pz)[k] < 0.0)
      (*pz)[k] += header1.BoxSize;
    while((*pz)[k] > header1.BoxSize)
      (*pz)[k] -= header1.BoxSize;
  }
  STOP;
  CHECK(sizeof(float)*3,*Np);

  //read velocities
  START;
  assert(fread(partChunk,(size_t) (*Np),3*sizeof(float),fp) == 3*sizeof(float));
  for(k=0;k<(*Np);++k) {
    (*vx)[k] = partChunk[k*3 + 0];
    (*vy)[k] = partChunk[k*3 + 1];
    (*vz)[k] = partChunk[k*3 + 2];
  }
  STOP;
  CHECK(sizeof(float)*3,*Np);
    
  //read IDs
  START;
  if(dummy_start == (*Np)*sizeof(int)) {
    assert(fread(idChunk_int,(size_t) (*Np),sizeof(int),fp) == sizeof(int));
    for(k=0;k<(*Np);++k)
      (*id)[k] = idChunk_int[k];    
    STOP;
    CHECK(sizeof(int),*Np);
  } else if(dummy_start == (*Np)*sizeof(long)) {
    assert(fread(idChunk_long,(size_t) (*Np),sizeof(long),fp) == sizeof(long));
    for(k=0;k<(*Np);++k)
      (*id)[k] = idChunk_long[k];    
    STOP;
    CHECK(sizeof(long),*Np);    
  } else {
    fprintf(stderr,"could not find correct size for IDs! %d for %ld particles\n",dummy_start,*Np);
    assert(0);
  }  
  
  //clean it all up
  fclose(fp);
  free(partChunk);
#undef START
#undef STOP
#undef CHECK
}

void write_parts_LGADGET(char fname_out[], char fname[], float *px, float *py, float *pz, float *vx, float *vy , float *vz, long *id, long Np)
{
  //vars and macros
  FILE *fp;
  struct io_header_1 header1;
  int dummy_start,dummy_stop,k;
  float *partChunk;
  long *idChunk;

  //get header from other file
#define RSTART fread(&dummy_start,sizeof(dummy_start),(size_t) 1,fp);
#define RSTOP fread(&dummy_stop,sizeof(dummy_stop),(size_t) 1,fp);
#define RCHECK(SZE,NP) assert(dummy_start == dummy_stop && (SZE)*(NP) == dummy_start);  
  fp = fopen_retry(fname,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname);
    assert(fp != NULL);
  }
  RSTART;
  fread(&header1,sizeof(header1),(size_t) 1,fp);
  RSTOP;
  RCHECK(256,1);  
  fclose(fp);
#undef RSTART
#undef RSTOP
#undef RCHECK
  
  //buffers and macros
  partChunk = (float*)malloc(sizeof(float)*Np*3);
  assert(partChunk != NULL);
  idChunk = (long*)partChunk;  
#define WSTART(SZE) dummy_start = (SZE); assert(1 == fwrite(&dummy_start,sizeof(dummy_start),1,fp));
#define WSTOP(SZE) dummy_stop = (SZE); assert(1 == fwrite(&dummy_stop,sizeof(dummy_stop),1,fp));
  
  //set header
  header1.num_files = 1;
  header1.npart[1] = Np;
  header1.npart[2] = 0;
  header1.npartTotal[1] = Np;
  header1.npartTotal[2] = 0;
  
  //write new file
  fp = fopen_retry(fname_out,"wb");
  if(fp == NULL) {
    fprintf(stderr,"could not open file '%s'!\n",fname_out);
    assert(fp != NULL);
  }
  
  //header
  WSTART(256);
  assert(1 == fwrite(&header1,sizeof(header1),1,fp));
  WSTOP(256);
  
  //positions                                                                                                                                                                                           
  for(k=0;k<Np;++k) {
    partChunk[k*3 + 0] = px[k];
    partChunk[k*3 + 1] = py[k];
    partChunk[k*3 + 2] = pz[k];
  }
  WSTART(3*sizeof(float)*Np);  
  assert(3*sizeof(float) == fwrite(partChunk,(size_t) Np,3*sizeof(float),fp));
  WSTOP(3*sizeof(float)*Np);  
    
  //velocities
  for(k=0;k<Np;++k) {
    partChunk[k*3 + 0] = vx[k];
    partChunk[k*3 + 1] = vy[k];
    partChunk[k*3 + 2] = vz[k];
  }
  WSTART(3*sizeof(float)*Np);  
  assert(3*sizeof(float) == fwrite(partChunk,(size_t) Np,3*sizeof(float),fp));
  WSTOP(3*sizeof(float)*Np);  
  
  //ids
  for(k=0;k<Np;++k)
    idChunk[k] = id[k];
  WSTART(sizeof(long)*Np);  
  assert(sizeof(long) == fwrite(idChunk,(size_t) Np,sizeof(long),fp));
  WSTOP(sizeof(long)*Np);  
    
  //clean it all up
  fclose(fp);
  free(partChunk);
}
