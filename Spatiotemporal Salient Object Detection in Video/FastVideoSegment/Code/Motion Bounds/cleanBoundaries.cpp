/** Function to remove boundaries that touch a given binary mask
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

/** Expected inputs:
			BoundaryMap: MxN logical matrix denoting bounds
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
			if( mxGetClassID( prhs[ 0 ] ) != mxLOGICAL_CLASS )
				mexErrMsgTxt( "Input matrix must be of type logical" );
			break;
		case 2:
			if( mxGetClassID( prhs[ 0 ] ) != mxLOGICAL_CLASS )
				mexErrMsgTxt( "Input matrix must be of type logical" );
			if( mxGetClassID( prhs[ 1 ] ) != mxLOGICAL_CLASS )
				mexErrMsgTxt( "Mask must be of type logical" );
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
	
	if( mxGetNumberOfDimensions( prhs[ 0 ] ) != 2 )
		mexErrMsgTxt( "Input matrix must be 2-dimensional" );
	if( ( nrhs == 2 ) && mxGetNumberOfDimensions( prhs[ 1 ] ) != 2 )
		mexErrMsgTxt( "Mask must be 2-dimensional" );

	bool *isBoundary = ( bool * )mxGetData( (mxArray *)prhs[ 0 ] );
	int height = mxGetM( (mxArray *)prhs[ 0 ] );
	int width = mxGetN( (mxArray *)prhs[ 0 ] );
	
	bool **remove = new bool*[ height ];
	for( int i = 0; i < height; i++ )
	{
		remove[ i ] = new bool[ width ];
		for( int j = 0; j < width; j++ )
			remove[ i ][ j ] = false;
	}
	
	/** Remove boundaries touching any of the edges */
	for( int i = 0; i < height; i++ )
	{
		if( isBoundary[ i ] )
			remove[ i ][ 0 ] = true;
		if( isBoundary[ (  width - 1 ) * height + i ] )
			remove[ i ][ width - 1 ] = true;
	}
	for( int j = 0; j < width; j++ )
	{
		if( isBoundary[ j * height ] )
			remove[ 0 ][ j ] = true;
		if( isBoundary[ j * height + height - 1 ] )
			remove[ height - 1 ][ j ] = true;
	}
	
	/** If a mask is provided, revove boundaries on the mask */
	if( nrhs == 2 )
	{
		bool *mask = ( bool * )mxGetData( (mxArray *)prhs[ 1 ] );
	
		int point;
		for( int i = 0; i < height; i++ )
		{
			for( int j = 0; j < width; j++ )
			{
				point = j * height + i;
				if( mask[ point ] && isBoundary[ point ] )
				{
					remove[ i ][ j ] = true;
				}
			}
		}
	}
	
	// Rightward scan
	for( int j = 1; j < width; j++ )
	{
		for( int i = 0; i < height; i++ )
		{
			if( isBoundary[ j * height + i ] )
			{	
				if( remove[ i ][ j - 1 ] )
					remove[ i ][ j ] = true;
				else
				{
					if( ( i > 0 ) && remove[ i - 1 ][ j - 1 ] )
						remove[ i ][ j ] = true;
					else if( ( i < height - 1 ) && remove[ i + 1 ][ j - 1 ] )
						remove[ i ][ j ] = true;
				}
			}
		}
	}
	
	// Leftward scan
	for( int j = width - 2; j > -1; j-- )
	{
		for( int i = 0; i < height; i++ )
		{
			if( isBoundary[ j * height + i ] )
			{
				if( remove[ i ][ j + 1 ] )
					remove[ i ][ j ] = true;
				else
				{
					if( ( i > 0 ) && remove[ i - 1 ][ j + 1 ] )
						remove[ i ][ j ] = true;
					else if( ( i < height - 1 ) && remove[ i + 1 ][ j + 1 ] )
						remove[ i ][ j ] = true;
				}
			}
		}
	}
	
	// Downward scan
	for( int i = 1; i < height; i++ )
	{
		for( int j = 0; j < width; j++ )
		{
			if( isBoundary[ j * height + i ] )
			{
				if( remove[ i - 1 ][ j ] )
					remove[ i ][ j ] = true;
				else
				{
					if( ( j > 0 ) && remove[ i - 1 ][ j - 1 ] )
						remove[ i ][ j ] = true;
					else if( ( j < width - 1 ) && remove[ i - 1 ][ j + 1 ] )
						remove[ i ][ j ] = true;
				}
			}
		}
	}
	
	// Upward scan
	for( int i = height - 2; i > -1; i-- )
	{
		for( int j = 0; j < width; j++ )
		{
			if( isBoundary[ j * height + i ] )
			{
				if( remove[ i + 1 ][ j ] )
					remove[ i ][ j ] = true;
				else
				{
					if( ( j > 0 ) && remove[ i + 1 ][ j - 1 ] )
						remove[ i ][ j ] = true;
					else if( ( j < width - 1 ) && remove[ i + 1 ][ j + 1 ] )
						remove[ i ][ j ] = true;
				}
			}
		}
	}
	
	mxArray *output = mxCreateLogicalMatrix( height, width );
	bool *outData = ( bool * )mxGetData( output );
	
	for( int i = 0; i < height; i++ )
	{
			for( int j = 0; j < width; j++ )
			{
				int point = j * height + i;
				outData[ point ] = isBoundary[ point ] & !remove[ i ][ j ];
				}
	}
	
	plhs[ 0 ] = output;
	
	for( int i = 0; i < height; i++ )
		delete [] remove[ i ];
	delete [] remove;

}

void printMatrix( const char* tag, unsigned char **matrix,
	short int height, short int width )
{
	mexPrintf( "%s\n", tag );
	for( short int i = 0; i < height; i++ )
	{
		for( short int j = 0; j < width; j++ )
			mexPrintf( "%i ", matrix[ i ][ j ] );
		mexPrintf( "\n" );
	}
	mexPrintf( "\n" );
}
