#include <matrix.h>
#include <mex.h>
#include <cmath> 
inline double min_of_three( double x, double y, double z )
{
	double min = y;

	if( x < min )
	{
		min = x;
	}
	if( z < min )
	{
		return z;
	}

	return min;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

    double *gradient;
    int M,N;
    int i,j;
    gradient = mxGetPr(prhs[0]);


    M = mxGetM(prhs[0]);//500行
    N =mxGetN(prhs[0]);//363列



//    mxArray *gf = mxCreateDoubleMatrix(M,N,mxREAL);
//    double *gfl = ( double * )mxGetData( gf );
	mxArray *energy = mxCreateDoubleMatrix(M,N,mxREAL);
    double *energyl = ( double * )mxGetData( energy );
//	for(j =0; j<N; j++)//列
//	{
//		for(i =0; i<M; i++)//行
//		{
//			if (j == 0)
//				gfl[j*M+i]=0;
//			else
//				gfl[j*M+i]=fabs( gradient[ j*M+i]-gradient[ (j-1)*M+i] );
//		}
//	}

	for(i =0; i<M; i++)//行
	{
			energyl[i]=gradient[ i];
	}

	for(j =1; j<N; j++)//列
	{
		for(i =0; i<M; i++)//行
		{
			if (i == 0)
			{
				energyl[j*M+i]= min_of_three( energyl[(j-1)*M+i],
											  energyl[(j-1)*M+i+1] ,
											  energyl[(j-1)*M+i+2]) +gradient[ j*M+i];
			}
			else if (i == 1)
			{
				double t =  min_of_three( energyl[(j-1)*M+i],
									      energyl[(j-1)*M+i+1],
										  energyl[(j-1)*M+i+2]);


				energyl[j*M+i]= min_of_three( t,
									      energyl[(j-1)*M+i-1],
										  energyl[(j-1)*M+i+2])+ gradient[ j*M+i];
			}
			else if (i == M-1)
			{
				energyl[j*M+i]= min_of_three( energyl[(j-1)*M+i],
											  energyl[(j-1)*M+i-1] ,
											  energyl[(j-1)*M+i-2]) + gradient[ j*M+i];
			}
			else if (i == M-2)
			{
				double t =  min_of_three( energyl[(j-1)*M+i],
										  energyl[(j-1)*M+i-1] ,
										  energyl[(j-1)*M+i-2]);

				energyl[j*M+i]= min_of_three( t,
									      energyl[(j-1)*M+i+1],
										  energyl[(j-1)*M+i-2])+ gradient[ j*M+i];
			}
			else
			{
				double t =  min_of_three( energyl[(j-1)*M+i],
										  energyl[(j-1)*M+i-1] ,
										  energyl[(j-1)*M+i+1]);

				energyl[j*M+i]= min_of_three( t,
									      energyl[(j-1)*M+i+2],
										  energyl[(j-1)*M+i-2])+ gradient[ j*M+i];
			}
		}
	}






    //for(i =0; i<M; i++)
    //{
    //     for(j =0; j<N; j++)
    //     {
    //         gfl[j*M+i]=gradient[ j*M+i];
    //     }
    // }
	//gfl[0] = M;
    plhs[ 0 ] = energy;
}

