#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <iostream>
#include <sys/time.h>

#define SIZE 16
float **makeArray2f(int width, int height);
void freeArray2f(float **arr, int width);
bool generateBMP(float **data, char *filePath, int dimension);
bool write2file(float **data,char *filePath, int dimension);



/*
    benchmark related vars.
*/

#define MR(mt,n,r,c,d)  mt->m[(n) * mt->mrows * mt->mcols * mt->mdeps + (r) * mt->mcols* mt->mdeps + (c) * mt->mdeps + (d)]

struct Mat {
  float* m;
  int mnums;
  int mrows;
  int mcols;
  int mdeps;
};

/* prototypes */
typedef struct Mat Matrix;

int newMat(Matrix* Mat, int mnums, int mrows, int mcols, int mdeps);
void clearMat(Matrix* Mat);
void mat_set(Matrix* Mat,int l,float z);
void mat_set_init(Matrix* Mat);
float jacobi(int n,Matrix* M1,Matrix* M2,Matrix* M3,Matrix* M4,Matrix* M5,Matrix* M6,Matrix* M7);
bool write2file3d_witharray(float *data,char *filePath,int row,int col , int depth);
float   omega=0.8;
Matrix  a,b,c,p,bnd,wrk1,wrk2;

/*
	vars related to pthread -as header 
*/
void *doStencil(void *);
void Barrier();
pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived=0,numWorkers = 0;       /* count of the number who have arrived */

// pthread structure para.
struct arguments {
    Matrix* a;
    Matrix* b;
    Matrix* c;
    Matrix* p;
    Matrix* bnd;
    Matrix* wrk1;
    Matrix* wrk2;
    int id;
    int nn;
};


/*
	pthread worker function (real part) -as header 
*/
void *doStencil(void *param) {
    struct arguments *args;
    Matrix* a;
    Matrix* b;
    Matrix* c;
    Matrix* p;
    Matrix* bnd;
    Matrix* wrk1;
    Matrix* wrk2;
    int id;
    int nn;
    
	// get argum.
    args = (struct arguments*)param;
    a = args->a;
    b = args->b;
    c = args->c;
    p = args->p;
    bnd = args->bnd;
    wrk1 = args->wrk1;
    wrk2 = args->wrk2;
    id = args->id;
    nn = args->nn;

	// devide data acc. i, j, k max 

    // impliment stencil 
    int    i,j,k,n,imax,jmax,kmax;
    float  gosa,s0,ss;

    imax= p->mrows-1;
    jmax= p->mcols-1;
    kmax= p->mdeps-1;

    for(n=0 ; n<nn ; n++){
        gosa = 0.0;

        for(i=1 ; i<imax; i++)
            for(j=1 ; j<jmax ; j++)
                for(k=1 ; k<kmax ; k++){
                    s0= MR(a,0,i,j,k)*MR(p,0,i+1,j,  k)
                        + MR(a,1,i,j,k)*MR(p,0,i,  j+1,k)
                        + MR(a,2,i,j,k)*MR(p,0,i,  j,  k+1)
                        + MR(b,0,i,j,k)
                        *( MR(p,0,i+1,j+1,k) - MR(p,0,i+1,j-1,k)
                        - MR(p,0,i-1,j+1,k) + MR(p,0,i-1,j-1,k) )
                        + MR(b,1,i,j,k)
                        *( MR(p,0,i,j+1,k+1) - MR(p,0,i,j-1,k+1)
                        - MR(p,0,i,j+1,k-1) + MR(p,0,i,j-1,k-1) )
                        + MR(b,2,i,j,k)
                        *( MR(p,0,i+1,j,k+1) - MR(p,0,i-1,j,k+1)
                        - MR(p,0,i+1,j,k-1) + MR(p,0,i-1,j,k-1) )
                        + MR(c,0,i,j,k) * MR(p,0,i-1,j,  k)
                        + MR(c,1,i,j,k) * MR(p,0,i,  j-1,k)
                        + MR(c,2,i,j,k) * MR(p,0,i,  j,  k-1)
                        + MR(wrk1,0,i,j,k);

                    ss= (s0*MR(a,3,i,j,k) - MR(p,0,i,j,k))*MR(bnd,0,i,j,k);

                    gosa+= ss*ss;
                    MR(wrk2,0,i,j,k)= MR(p,0,i,j,k) + omega*ss;
                }

        for(i=1 ; i<imax ; i++)
            for(j=1 ; j<jmax ; j++)
                for(k=1 ; k<kmax ; k++)
                    MR(p,0,i,j,k)= MR(wrk2,0,i,j,k);
        
    } /* end n loop */
	return NULL;
}

/*
	Main Function 
*/
int main(int argc, char** argv) {
    std::cin>> numWorkers;
	
    int    nn;
    int    imax,jmax,kmax,mimax,mjmax,mkmax,msize[3];
    float  gosa;

    //   scanf("%d", &msize[0]);
    //   scanf("%d", &msize[1]);
    //   scanf("%d", &msize[2]);
    //   scanf("%d", &nn);
    
    mimax= 16;
    mjmax= 16;
    mkmax= 128;
    nn=10;
    imax= mimax-1;
    jmax= mjmax-1;
    kmax= mkmax-1;

    /*
    *    Initializing matrixes
    */
    newMat(&p,1,mimax,mjmax,mkmax);
    newMat(&bnd,1,mimax,mjmax,mkmax);
    newMat(&wrk1,1,mimax,mjmax,mkmax);
    newMat(&wrk2,1,mimax,mjmax,mkmax);
    newMat(&a,4,mimax,mjmax,mkmax);
    newMat(&b,3,mimax,mjmax,mkmax);
    newMat(&c,3,mimax,mjmax,mkmax);

    mat_set_init(&p);
    mat_set(&bnd,0,1.0);
    mat_set(&wrk1,0,0.0);
    mat_set(&wrk2,0,0.0);
    mat_set(&a,0,1.0);
    mat_set(&a,1,1.0);
    mat_set(&a,2,1.0);
    mat_set(&a,3,1.0/6.0);
    mat_set(&b,0,0.0);
    mat_set(&b,1,0.0);
    mat_set(&b,2,0.0);
    mat_set(&c,0,1.0);
    mat_set(&c,1,1.0);
    mat_set(&c,2,1.0);

    /*
    *    do sequential
    */
    gosa = jacobi(nn,&a,&b,&c,&p,&bnd,&wrk1,&wrk2);

    Matrix* pp=&p;
    write2file3d_witharray(pp->m,"result_seq.txt",pp->mrows,pp->mcols,pp->mdeps);

    /*
        after seq. initialize
    */
    mat_set_init(&p);
    mat_set(&bnd,0,1.0);
    mat_set(&wrk1,0,0.0);
    mat_set(&wrk2,0,0.0);
    mat_set(&a,0,1.0);
    mat_set(&a,1,1.0);
    mat_set(&a,2,1.0);
    mat_set(&a,3,1.0/6.0);
    mat_set(&b,0,0.0);
    mat_set(&b,1,0.0);
    mat_set(&b,2,0.0);
    mat_set(&c,0,1.0);
    mat_set(&c,1,1.0);
    mat_set(&c,2,1.0);

    /*
        do parallel computing 
    */
    pthread_t workerid[numWorkers];
    struct arguments arg[numWorkers];
    // create threads
    for(int i = 0; i < numWorkers; i++) {
        // set args
        arg[i].nn = nn;
        arg[i].id = i;
        arg[i].a = &a;
        arg[i].b = &b;
        arg[i].c = &c;
        arg[i].p = &p;
        arg[i].bnd = &bnd;
        arg[i].wrk1 = &wrk1;
        arg[i].wrk2 = &wrk2;
		std::cout<< "starting thread "<<i<<std::endl;
        pthread_create(&workerid[i], NULL, doStencil, (void*)&arg[i]);
    }
    // Wait for all threads to complete
    for(int i = 0; i < numWorkers; i++) {
        pthread_join(workerid[i], NULL);
    }
    printf("Done\n");

    pp=&p;
    write2file3d_witharray(pp->m,"result_parallel.txt",pp->mrows,pp->mcols,pp->mdeps);

    //printf("%.6f\n",gosa);

    /*
    *   Matrix free
    */ 
    clearMat(&p);
    clearMat(&bnd);
    clearMat(&wrk1);
    clearMat(&wrk2);
    clearMat(&a);
    clearMat(&b);
    clearMat(&c);

	return 0;
}



/*
 * dynamic 2D array of floats allocation 
 */
float **makeArray2f(int width, int height) {
    size_t rowSize = sizeof(float*) * width;
    size_t colSize = sizeof(float) * height;
    float **result = (float**)malloc(rowSize);
    for(int i = 0; i < width; i++) {
        result[i] = (float*)malloc(colSize);
    }
    return result;
}

/*
 * Frees memory from a 2D float array
 */
void freeArray2f(float **arr, int width) {
    for(int i = 0; i < width; i++) {
        free(arr[i]);
    }
    free(arr);
}

/*
 * Writes data to a 24-bit BMP file.
 * assumes input is dim x dim array
 */
bool generateBMP(float **data, char *filePath, int dimension) {
    // Open file.  Overwrites existing file if it exists
    FILE *bitFile;
    if(!(bitFile = fopen(filePath, "wb"))) return false;
    // Write header
    int rowSize = 4 * ((int)((24 * dimension + 31) / 32));
    int fileSize = 54 + rowSize * dimension;
    char hd1[2] = {0x42, 0x4D};     // Magic Numbers
    short hd2[2] = {0x0, 0x0};      // Application Specific Values
    int offset = 54;
    int dibHeaderData1[3] = {40, dimension, dimension}; // header size, width, height
    short dibHeaderData2[2] = {1, 24};  // color planes, bits/pixel
    int dibHeaderData3[6] = {0, fileSize - offset, 2835, 2835, 0, 0};
                            // RGB, size of data, horz. resolution, vert res
    fwrite(hd1, 1, 2, bitFile);
    fwrite(&fileSize, 4, 1, bitFile);
    fwrite(hd2, 2, 2, bitFile);
    fwrite(&offset, 4, 1, bitFile);
    fwrite(dibHeaderData1, 4, 3, bitFile);
    fwrite(dibHeaderData2, 2, 2, bitFile);
    fwrite(dibHeaderData3, 4, 6, bitFile);
    // Write pixel data, starts from bottom row
    // each row should begin at a location in file that is a multiple of 4
    for(int i = dimension - 1; i >= 0; i--) {
        for(int j = 0; j < dimension; j++) {
            // BGR order
            char color[3] = {(char)(255 * (1 - data[i][j])),
                            0,
                            (char)(255 * data[i][j])};
            fwrite(color, 1, 3, bitFile);
        }
        // Check that next row will be at multiple of 4
        int remainder = (dimension * 3) % 4;
        if(remainder > 0) {
            char pad = 0;
            fwrite(&pad, remainder, 1, bitFile);
        }
    }
    fclose(bitFile);
    return true;
}
bool write2file(float **data,char *filePath, int dimension){

	FILE * results = fopen(filePath, "w");
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			fprintf(results, "%.6f ", data[i][j]);
			}
				fprintf(results, "\n");
			}
	return true ;
}
void Barrier() {
	pthread_mutex_lock(&barrier);
	numArrived++;
	if (numArrived == numWorkers) {
		//std::cout<< "go!"<<std::endl;
		numArrived = 0;
		pthread_cond_broadcast(&go);
	} else
		pthread_cond_wait(&go, &barrier);
	pthread_mutex_unlock(&barrier);
}

int
newMat(Matrix* Mat, int mnums,int mrows, int mcols, int mdeps)
{
  Mat->mnums= mnums;
  Mat->mrows= mrows;
  Mat->mcols= mcols;
  Mat->mdeps= mdeps;
  Mat->m= NULL;
  Mat->m= (float*) 
    malloc(mnums * mrows * mcols * mdeps * sizeof(float));
  
  return(Mat->m != NULL) ? 1:0;
}

void
clearMat(Matrix* Mat)
{
  if(Mat->m)
    free(Mat->m);
  Mat->m= NULL;
  Mat->mnums= 0;
  Mat->mcols= 0;
  Mat->mrows= 0;
  Mat->mdeps= 0;
}

void
mat_set(Matrix* Mat, int l, float val)
{
  int i,j,k;

    for(i=0; i<Mat->mrows; i++)
      for(j=0; j<Mat->mcols; j++)
        for(k=0; k<Mat->mdeps; k++)
          MR(Mat,l,i,j,k)= val;
}

void
mat_set_init(Matrix* Mat)
{
  int  i,j,k;

  for(i=0; i<Mat->mrows; i++)
    for(j=0; j<Mat->mcols; j++)
      for(k=0; k<Mat->mdeps; k++)
        MR(Mat,0,i,j,k)= (float)(i*i)
          /(float)((Mat->mrows - 1)*(Mat->mrows - 1));
}

float
jacobi(int nn, Matrix* a,Matrix* b,Matrix* c,
       Matrix* p,Matrix* bnd,Matrix* wrk1,Matrix* wrk2)
{
  int    i,j,k,n,imax,jmax,kmax;
  float  gosa,s0,ss;

  imax= p->mrows-1;
  jmax= p->mcols-1;
  kmax= p->mdeps-1;

  for(n=0 ; n<nn ; n++){
    gosa = 0.0;

    for(i=1 ; i<imax; i++)
      for(j=1 ; j<jmax ; j++)
        for(k=1 ; k<kmax ; k++){
          s0= MR(a,0,i,j,k)*MR(p,0,i+1,j,  k)
            + MR(a,1,i,j,k)*MR(p,0,i,  j+1,k)
            + MR(a,2,i,j,k)*MR(p,0,i,  j,  k+1)
            + MR(b,0,i,j,k)
             *( MR(p,0,i+1,j+1,k) - MR(p,0,i+1,j-1,k)
              - MR(p,0,i-1,j+1,k) + MR(p,0,i-1,j-1,k) )
            + MR(b,1,i,j,k)
             *( MR(p,0,i,j+1,k+1) - MR(p,0,i,j-1,k+1)
              - MR(p,0,i,j+1,k-1) + MR(p,0,i,j-1,k-1) )
            + MR(b,2,i,j,k)
             *( MR(p,0,i+1,j,k+1) - MR(p,0,i-1,j,k+1)
              - MR(p,0,i+1,j,k-1) + MR(p,0,i-1,j,k-1) )
            + MR(c,0,i,j,k) * MR(p,0,i-1,j,  k)
            + MR(c,1,i,j,k) * MR(p,0,i,  j-1,k)
            + MR(c,2,i,j,k) * MR(p,0,i,  j,  k-1)
            + MR(wrk1,0,i,j,k);

          ss= (s0*MR(a,3,i,j,k) - MR(p,0,i,j,k))*MR(bnd,0,i,j,k);

          gosa+= ss*ss;
          MR(wrk2,0,i,j,k)= MR(p,0,i,j,k) + omega*ss;
        }

    for(i=1 ; i<imax ; i++)
      for(j=1 ; j<jmax ; j++)
        for(k=1 ; k<kmax ; k++)
          MR(p,0,i,j,k)= MR(wrk2,0,i,j,k);
    
  } /* end n loop */

  return(gosa);
}

bool write2file3d_witharray(float *data,char *filePath,int row,int col ,int depth){

	FILE * results = fopen(filePath, "w");
    for (int k = 0; k < depth; k++) {
        fprintf(results, "k= %d\n" , k);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                fprintf(results, "%.6f ", data[k*row*col+i*col+j]);
            }
            fprintf(results, "\n");
        }
        
        fprintf(results, "\n");
    }
	return true ;
}


