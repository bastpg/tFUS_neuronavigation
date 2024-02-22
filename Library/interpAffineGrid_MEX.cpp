


#include "mex.h"
#include "matrix.h"
#include <math.h>


#define DEBUG  0



void interpAffineGrid( double* x , double* y , double* z , double* DATA , double* x2 , double* y2 , double* z2 , double* R , 
                        int nx , int ny , int nz , int nx2 , int ny2 , int nz2 , double *DATA2 )
{

    double Lx = x[nx-1] - x[0];
    double Ly = y[ny-1] - y[0];
    double Lz = z[nz-1] - z[0];

    double tPOS[4] , iPOS[4];    
    
    for( int ii = 0 ; ii < nx2 ; ii++ )
    {
        
        if( DEBUG == 1 )
            printf( "  ii = %d / %d\n" , ii , nx2 );
    
        for( int jj = 0 ; jj < ny2 ; jj++ )
        {

            for( int kk = 0 ; kk < nz2 ; kk++ )
            {
                
                // compute affine transformed grid
                iPOS[0] = x2[ii];
                iPOS[1] = y2[jj];
                iPOS[2] = z2[kk];
                iPOS[3] = 1.0;
                
                for( int ll = 0 ; ll < 4  ; ll++ ){
                    tPOS[ll] = 0.0;
                    for( int mm = 0 ; mm < 4  ; mm++ ){
                        tPOS[ ll ]  +=   R[ ll + mm * 4 ] * iPOS[mm];
                    }                
                }
                
                
                if( DEBUG == 1 ){
                    if( ii == 0 & jj == 0 & kk == 0 ){
                        printf( "*** ii = %d   jj = %d  kk = %d  ***\n" , ii , jj , kk );
                        printf( "   iPOS[0] = %f     tPOS[0] = %f\n" , iPOS[0] , tPOS[0] );
                        printf( "   iPOS[1] = %f     tPOS[1] = %f\n" , iPOS[1] , tPOS[1] );
                        printf( "   iPOS[2] = %f     tPOS[2] = %f\n" , iPOS[2] , tPOS[2] );
                        printf( "   iPOS[3] = %f     tPOS[3] = %f\n\n" , iPOS[3] , tPOS[3] );
                    }
                }                    
                
                double ii2 = ( tPOS[0] - x[0] ) / Lx * ( nx - 1 );    // fractional indices in the orginal grid (= Matlab indices - 1)
                double jj2 = ( tPOS[1] - y[0] ) / Ly * ( ny - 1 );
                double kk2 = ( tPOS[2] - z[0] ) / Lz * ( nz - 1 );

                int ii2_ = int(ii2);
                int jj2_ = int(jj2);
                int kk2_ = int(kk2);

                if( ii2_ >= 0  &  ii2_ < nx-1  &  jj2_ >= 0  &  jj2_ < ny-1  &  kk2_ >= 0  &  kk2_ < nz-1 )
                {

                    double dx = ii2 -  double( ii2_ );
                    double dy = jj2 -  double( jj2_ );
                    double dz = kk2 -  double( kk2_ );
                    
                    DATA2[ ii + jj*nx2 + kk*nx2*ny2 ] = DATA[ ii2_   +     jj2_*nx  +      kk2_*nx*ny ] * (1 - dx) * (1 - dy) * (1 - dz)   +    
                                                        DATA[ ii2_+1 +     jj2_*nx  +      kk2_*nx*ny ] * (  dx  ) * (1 - dy) * (1 - dz)   +   
                                                        DATA[ ii2_   + (jj2_+1)*nx  +      kk2_*nx*ny ] * (1 - dx) * (  dy  ) * (1 - dz)   + 
                                                        DATA[ ii2_   +     jj2_*nx  +  (kk2_+1)*nx*ny ] * (1 - dx) * (1 - dy) * (  dz  )   +  
                                                        DATA[ ii2_+1 + (jj2_+1)*nx  +      kk2_*nx*ny ] * (  dx  ) * (  dy  ) * (1 - dz)   + 
                                                        DATA[ ii2_+1 +     jj2_*nx  +  (kk2_+1)*nx*ny ] * (  dx  ) * (1 - dy) * (  dz  )   +  
                                                        DATA[ ii2_   + (jj2_+1)*nx  +  (kk2_+1)*nx*ny ] * (1 - dx) * (  dy  ) * (  dz  )   + 
                                                        DATA[ ii2_+1 + (jj2_+1)*nx  +  (kk2_+1)*nx*ny ] * (  dx  ) * (  dy  ) * (  dz  );


                }  // if
                         
                
            }  // for x
        }  // for y
    }  // for z
    
}






















void mexFunction( int nlhs , mxArray *plhs[] , int nrhs , const mxArray *prhs[] )
{
    if( nrhs != 8 ){
        printf("Usage:\n\tDATA2 = interpAffineGrid( x , y , z , DATA , x2 , y2 , z2 , R )\n");                
        mexErrMsgTxt("Wrong number of input parameters.");        
    }

    // STRUCTURES
    double *x = (double *)mxGetPr(prhs[0]);  
    double *y = (double *)mxGetPr(prhs[1]);  
    double *z = (double *)mxGetPr(prhs[2]);  
    
    double *DATA = (double *)mxGetPr(prhs[3]);  

    double *x2 = (double *)mxGetPr(prhs[4]);  
    double *y2 = (double *)mxGetPr(prhs[5]);  
    double *z2 = (double *)mxGetPr(prhs[6]);  
    
    double *R = (double *)mxGetPr(prhs[7]);  

    int nx = (int)mxGetNumberOfElements(prhs[0]);
    int ny = (int)mxGetNumberOfElements(prhs[1]);
    int nz = (int)mxGetNumberOfElements(prhs[2]);
    
    int nx2 = (int)mxGetNumberOfElements(prhs[4]);
    int ny2 = (int)mxGetNumberOfElements(prhs[5]);
    int nz2 = (int)mxGetNumberOfElements(prhs[6]);
            
    if( DEBUG == 1 ){
        printf( "  nx = %d  ny = %d  nz = %d  nx2 = %d  ny2 = %d  nz2 = %d\n\n" , nx , ny , nz , nx2 , ny2 , nz2 );
         printf( "\n   R = [  %f  %f  %f  %f\n" , R[0+0*4] , R[0+1*4] , R[0+2*4] , R[0+3*4] );   
         printf( "          %f  %f  %f  %f\n"   , R[1+0*4] , R[1+1*4] , R[1+2*4] , R[1+3*4] );   
         printf( "          %f  %f  %f  %f\n"   , R[2+0*4] , R[2+1*4] , R[2+2*4] , R[2+3*4] );    
         printf( "          %f  %f  %f  %f  ]\n\n" , R[3+0*4] , R[3+1*4] , R[3+2*4] , R[3+3*4] ); 
    }
    
    
    // CONSISTENCY CHECK
    if( (nx*ny*nz) != (int)mxGetNumberOfElements(prhs[3])  |  (int)mxGetNumberOfElements(prhs[7]) != 16  ){
        mexErrMsgTxt("Inconsistent dimension of input structures." );
    }
    
    // ALLOCATE MEMORY
    double *DATA2;
    
    plhs[0] = mxCreateDoubleMatrix( nx2 * ny2 * nz2 , 1 , mxREAL );
    DATA2 = (double *)mxGetPr(plhs[0]);
    
    if(nlhs>2){
        mexErrMsgTxt("Too many output variables.");
    }
    
    
    // PERFORM COMPUTATION        
    interpAffineGrid( x , y , z , DATA , x2 , y2 , z2 , R ,
            nx , ny , nz , nx2 , ny2 , nz2 , DATA2 );

}



























