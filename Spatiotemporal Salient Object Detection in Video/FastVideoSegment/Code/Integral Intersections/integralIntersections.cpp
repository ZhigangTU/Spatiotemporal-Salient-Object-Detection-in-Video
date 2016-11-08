/** Function that computes the integral intersections algorithm
%
%    Copyright (C) 2013  Anestis Papazoglou
%
%    Copyright (C) 2013  Anestis Papazoglou
%
%    You can redistribute and/or modify this software for non-commercial use
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    For commercial use, contact the author for licensing options.
%
%    Contact: a.papazoglou@sms.ed.ac.uk

/** TODO: Function is not fully optimised.
	Modulous checks can be replaced by flip-flop like operations.
	Furthermore, function can be threaded for speedup.
*/

#include <matrix.h>
#include <mex.h>

#define min( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )

void printMatrix( const char* tag, unsigned char **matrix,
	short int height, short int width );

void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[] )
{
	switch( nrhs )
	{
		case 1:
			if( mxGetClassID( prhs[ 0 ] ) != mxUINT8_CLASS )
				mexErrMsgTxt( "Input matrix must be of type uint8" );
			break;
		default:
			mexErrMsgTxt( "Invalid number of input arguments" );
	}

	switch( nlhs )
	{
		case 0:
			break;
		case 1:
			break;
		default:
			mexErrMsgTxt( "Invalid number of output arguments" );
	}
	
	if( mxGetNumberOfDimensions( prhs[ 0 ] ) > 2 )
		mexErrMsgTxt( "Input matrix must be 2-dimensional" );
	
	mxArray *input = (mxArray *)prhs[ 0 ];
	short int height = mxGetM( input );
	short int width = mxGetN( input );
	mxArray *output =
		mxCreateNumericMatrix( height, width, mxUINT8_CLASS, mxREAL );

	unsigned char *outputData = ( unsigned char * )mxGetData( output );
	plhs[ 0 ] = output;

	unsigned char **integralMatrix;

	integralMatrix = new unsigned char*[ height ];
	for( short int i = 0; i < height; i++ )
		integralMatrix[ i ] = new unsigned char[ width ];

	unsigned char integral;

	/** Horizontal integral line **/
	for( short int i = 0; i < height; i++ )
	{
		integralMatrix[ i ][ 0 ] = 0;
		for( short int j = 1; j < width; j++ )
		{
			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ ( j - 1 ) * height + i ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i ][ j - 1 ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i ][ j - 1 ];
		}
	}

	/** Horizontal lines check **/
	for( int i = 0; i < height; i++ )
	{
		for( int j = 1; j < width; j++ )
		{
			if( integralMatrix[ i ][ j - 1 ] % 2 )
				outputData[ j * height + i ]++;
		}
	}

	for( int i = 0; i < height; i++ )
	{
		for( int j = 0; j < width - 1; j++ )
		{
			if( ( integralMatrix[ i ][ width - 1 ] - integralMatrix[ i ][ j ] ) % 2 )
				outputData[ j * height + i ]++;
		}
	}

	/** Vertical integral line **/
	for( short int j = 0; j < width; j++ )
	{
		integralMatrix[ 0 ][ j ] = 0;
		for( short int i = 1; i < height; i++ )
		{
			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + ( i - 1 ) ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j ];
		}
	}

	/** Vertical lines check **/
	for( int i = 1; i < height; i++ )
	{
		for( int j = 0; j < width; j++ )
		{
			if( integralMatrix[ i - 1 ][ j ] % 2 )
				outputData[ j * height + i ]++;
		}
	}

	for( int i = 0; i < height - 1; i++ )
	{
		for( int j = 0; j < width; j++ )
		{
			if( ( integralMatrix[ height - 1 ][ j ] - integralMatrix[ i ][ j ] ) % 2 )
				outputData[ j * height + i ]++;
		}
	}

	/** Diagonal downward integral line **/
	for( short int start = 0; start < height; start++ )
	{
		short int stop = min( height - start, width );		
		short int i = start;
		short int j = 0;
		
		integralMatrix[ i ][ j ] = 0;
		for( short int step = 1; step < stop; step++ )
		{
			i = start + step;
			j = step;

			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ ( j - 1 ) * height + ( i - 1 ) ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j - 1 ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j - 1 ];
		}
	}

	for( short int start = 1; start < width; start++ )
	{
		short int stop = min( height, width - start );
		short int i = 0;
		short int j = start;

		integral = ( ( unsigned char* )( mxGetData( input ) ) )[ j * height ];
		integralMatrix[ i ][ j ] = 0;
		for( short int step = 1; step < stop; step++ )
		{
			i = step;
			j = start + step;

			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ ( j - 1 ) * height + ( i - 1 ) ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j - 1 ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i - 1 ][ j - 1 ];
		}
	}

	/** Diagonal downward lines check **/
	for( int i = 1; i < height; i++ )
	{
		for( int j = 1; j < width; j++ )
		{
			if( integralMatrix[ i - 1 ][ j - 1 ] % 2 )
				outputData[ j * height + i ]++;
		}
	}

	for( int i = 0; i < height - 1; i++ )
	{
		for( int j = 0; j < width - 1; j++ )
		{
			short int endx;
			short int endy;

			if( height - i < width - j )
			{
				endx = height - 1;
				endy = j + height - 1 - i;
			}
			else
			{
				endx = i + width - 1 - j;
				endy = width - 1;
			}

			if( ( integralMatrix[ endx ][ endy ] - integralMatrix[ i ][ j ] ) % 2 )
				outputData[ j * height + i ]++;
		}
	}

	/** Diagonal upward integral line **/
	for( short int start = height - 1; start > -1; start-- )
	{
		short int stop = min( start + 1, width );
		short int i = start;
		short int j = 0;

		integralMatrix[ i ][ j ] = 0;
		for( short int step = 1; step < stop; step++ )
		{
			i = start - step;
			j = step;

			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ ( j - 1 ) * height + ( i + 1 ) ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i + 1 ][ j - 1 ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i + 1 ][ j - 1 ];
		}
	}

	for( short int start = 1; start < width; start++ )
	{
		short int stop = min( height, width - start ); 		
		short int i = height - 1;
		short int j = start;

		integralMatrix[ i ][ j ] = 0;
		for( short int step = 1; step < stop; step++ )
		{
			i = height - 1 - step;
			j = start + step;

			if( ( ( unsigned char* )( mxGetData( input ) ) )[ j * height + i ] > ( ( unsigned char* )( mxGetData( input ) ) )[ ( j - 1 ) * height + ( i + 1 ) ] )
				integralMatrix[ i ][ j ] = integralMatrix[ i + 1 ][ j - 1 ] + 1;
			else
				integralMatrix[ i ][ j ] = integralMatrix[ i + 1 ][ j - 1 ];
		}
	}

	/** Diagonal upward lines check **/
	for( int i = height - 2; i > - 1; i-- )
	{
		for( int j = 1; j < width; j++ )
		{
			if( integralMatrix[ i + 1 ][ j - 1 ] % 2 )
				outputData[ j * height + i ]++;
		}
	}
	
	for( int i = height - 1; i > 0; i-- )
	{
		for( int j = 0; j < width - 1; j++ )
		{
			short int endx;
			short int endy;

			if( i < width - 1 - j )
			{
				endx = 0;
				endy = j + i;
			}
			else
			{
				endx = i - width + 1 + j;
				endy = width - 1;
			}

			if( ( integralMatrix[ endx ][ endy ] - integralMatrix[ i ][ j ] ) % 2 )
				outputData[ j * height + i ]++;
		}
	}

	/** Deallocate memory */
	for( short int i = 0; i < height; i++ )
	{
		delete [] integralMatrix[ i ];
	}
	delete [] integralMatrix;

}

void printMatrix( const char* tag, unsigned char **matrix,
	short int height, short int width )
{
	mexPrintf( "%s\n[ ", tag );
	for( short int i = 0; i < height; i++ )
	{
		for( short int j = 0; j < width; j++ )
			mexPrintf( "%i ", matrix[ i ][ j ] );
		mexPrintf( ",\n" );
	}
	mexPrintf( " ]\n" );
}
