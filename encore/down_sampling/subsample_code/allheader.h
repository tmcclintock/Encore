#include <stdio.h>

#ifndef _SUBSAMP_
#define _SUBSAMP_

#define MAX_FILENAME 2048

extern char SIMTYPE[1024];

/* in read_LGADGET.c */
long read_num_chunks_LGADGET(char fname[]);
long read_num_parts_LGADGET(char fname[]);
float read_period_length_LGADGET(char fname[]);
void read_parts_LGADGET(char fname[], float **px, float **py, float **pz, float **vx, float **vy , float **vz, long **id, long *Np);
void write_parts_LGADGET(char fname_out[], char fname[], float *px, float *py, float *pz, float *vx, float *vy , float *vz, long *id, long Np);

/* in io.c */
long read_num_chunks_file(char fname[]);
void read_parts_file(char fname[], float **px, float **py, float **pz, float **vx, float **vy , float **vz, long **id, long *Np);
void write_parts_file(char fname_out[], char fname[], float *px, float *py, float *pz, float *vx, float *vy , float *vz, long *id, long Np);
float read_period_length_file(char fname[]);
long read_num_parts_file(char fname[]);
FILE *fopen_retry(const char *filename, const char *mode);

#endif /* _SUBSAMP_ */
