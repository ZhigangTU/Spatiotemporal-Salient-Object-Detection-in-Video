/** Function to return superpixel labellings given superpixel boundaries
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
%    Contact: a.papazoglou@sms.ed.ac.uk */

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
			if( mxGetClassID( prhs[ 0 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Input matrix must be of type uint16" );
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

	int height = mxGetM( (mxArray *)prhs[ 0 ] );
	int width = mxGetN( (mxArray *)prhs[ 0 ] );
	unsigned short *inData = ( unsigned short * )mxGetData( prhs[ 0 ] );
	
	mxArray *output = mxDuplicateArray( prhs[0] );
	unsigned short *outData = ( unsigned short * )mxGetData( output );
	
	for( int i = 0; i < height; i++ )
	{
			for( int j = 0; j < width; j++ )
			{
				int point = j * height + i;
				if( ~outData[ point ] )
				{
					/** Check above and below */
					if( ( i > 0 ) && inData[ j * height + i - 1 ] )
						outData[ point ] = outData[ j * height + i - 1 ];
					else if( ( i < height - 1 ) && outData[ j * height + i + 1 ] )
						outData[ point ] = inData[ j * height + i + 1 ];
					/** Check right and left */
					else if( ( j > 0 ) && inData[ ( j - 1 ) * height + i ] )
						outData[ point ] = outData[ ( j - 1 ) * height + i ];
					else if( ( j < width - 1 ) && outData[ ( j + 1 ) * height + i ] )
						outData[ point ] = inData[ ( j + 1 ) * height + i ];
					/** Check diagonally */
					else if( ( i > 0 ) && ( j > 0 ) && inData[ ( j - 1 ) * height + i - 1 ] )
						outData[ point ] = outData[ ( j - 1 ) * height + i - 1 ];
					else if( ( i > 0 ) && ( j < width - 1 ) && inData[ ( j + 1 ) * height + i - 1 ] )
						outData[ point ] = outData[ ( j + 1 ) * height + i - 1 ];
					else if( ( i < height - 1 ) && ( j > 0 ) && inData[ ( j - 1 ) * height + i + 1 ] )
						outData[ point ] = outData[ ( j - 1 ) * height + i + 1 ];
					else if( ( i < height - 1 ) && ( j < width - 1 ) && inData[ ( j + 1 ) * height + i + 1 ] )
						outData[ point ] = outData[ ( j + 1 ) * height + i + 1 ];
					else
						mexPrintf( "Warning: point (%i,%i) could not be filled in - all neighbours are unlabeled", i, j );
				}
			}
	}
	
	plhs[ 0 ] = output;

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
