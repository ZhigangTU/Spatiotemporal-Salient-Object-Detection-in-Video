#ifndef CTENSORMATLAB_H
#define CTENSORMATLAB_H

#include <iostream>

#include <matrix.h>
#include <mex.h>

#include "CTensor.h"

template <class T> class CTensorMatlab: public CTensor< T >
{
	public:
		// standard constructor
		inline CTensorMatlab(): CTensor< T >() {}
		// constructor
		inline CTensorMatlab( const int xSize, const int ySize, const int zSize ): CTensor< T >( xSize, ySize, zSize ) {}
		// copy constructor
		CTensorMatlab( const CTensor<T> &aCopyFrom ): CTensor< T >( aCopyFrom ) {}
		// constructor with implicit filling
		CTensorMatlab( const int xSize, const int ySize, const int zSize, const T value ): CTensor< T >( xSize, ySize, zSize, value ) {}

		void copyFromMxArrayImage( const mxArray *input )
		{
			int dimensions = mxGetNumberOfDimensions( input );
			const int *size = mxGetDimensions( input );
			int xSize, ySize, zSize;

			if( dimensions == 2 )
			{
				xSize = size[ 0 ];
				ySize = size[ 1 ];
				zSize = 1;
			}
			else if( dimensions == 3 )
			{
				xSize = size[ 0 ];
				ySize = size[ 1 ];
				zSize = size[ 2 ];

				if( zSize != 3 )
				{
					//throw MxMatrixSizeIncompatible( xSize, ySize, zSize );
					std::cout << "crap" << std::endl;
					return;
				}
			}
			else
			{
				//throw MxMatrixSizeIncompatible( dimensions );
				std::cout << "crap" << std::endl;
				return;
			}

			if( xSize != CTensor< T >::mXSize || ySize != CTensor< T >::mYSize || zSize != CTensor< T >::mZSize )
			{
				delete [] CTensor< T >::mData;

				CTensor< T >::mXSize = xSize;
				CTensor< T >::mYSize = ySize;
				CTensor< T >::mZSize = zSize;
				CTensor< T >::mData = new T[ CTensor< T >::mXSize * CTensor< T >::mYSize * CTensor< T >::mZSize ];
			}

			// We assume the mxArray contains UINT8
			unsigned char *inputData = ( unsigned char * )mxGetData( input );
			int inputSize = mxGetNumberOfElements( input );
			for( int i = 0; i < inputSize; i++ )
			{
				CTensor< T >::mData[ i ] = ( T )inputData[ i ];
			}
		}

		void copyFromMxArrayOpticalFlow( mxArray *input )
		{
			int dimensions = mxGetNumberOfDimensions( input );
			const int *size = mxGetDimensions( input );

			int xSize, ySize, zSize;
			if( dimensions == 3 )
			{
				xSize = size[ 0 ];
				ySize = size[ 1 ];
				zSize = size[ 2 ];

				if( zSize != 2 )
				{
					//throw MxMatrixSizeIncompatible( xSize, ySize, zSize );
					std::cout << "crap" << std::endl;
					return;
				}
			}
			else
			{
				//throw MxMatrixSizeIncompatible( dimensions );
				std::cout << "crap" << std::endl;
				return;
			}

			if( xSize != CTensor< T >::mXSize || ySize != CTensor< T >::mYSize || zSize != CTensor< T >::mZSize )
			{
				delete [] CTensor< T >::mData;

				CTensor< T >::mXSize = xSize;
				CTensor< T >::mYSize = ySize;
				CTensor< T >::mZSize = zSize;
				CTensor< T >::mData = new T[ CTensor< T >::mXSize * CTensor< T >::mYSize * CTensor< T >::mZSize ];
			}

			// We assume the mxArray contains INT16
			short *inputData = ( short * )mxGetData( input );
			int inputSize = mxGetNumberOfElements( input );
			for( int i = 0; i < inputSize; i++ )
			{
				CTensor< T >::mData[ i ] = ( T )inputData[ i ];
			}
		}

		void copyToMxArray( mxArray *output )
		{
			T *data = ( T * )mxGetData( output );
			
			unsigned int size = CTensor< T >::mXSize * CTensor< T >::mYSize * CTensor< T >::mZSize;
			
			for( unsigned int i; i < size; i++ )
			{
				data[ i ] = CTensor< T >::mData[ i ];
			}
			
		}

};

struct MxMatrixSizeIncompatible
{
	MxMatrixSizeIncompatible( int dimensions )
	{
    std::cerr << "Exception MxMatrixSizeIncompatible: dimensions = " << dimensions << " (must be 2 or 3)" << std::endl;
  }
  MxMatrixSizeIncompatible( int x, int y, int z )
	{
    std::cerr << "Exception MxMatrixSizeIncompatible: " << x << "x" << y << "x" << z << std::endl;
  }
};

#endif
