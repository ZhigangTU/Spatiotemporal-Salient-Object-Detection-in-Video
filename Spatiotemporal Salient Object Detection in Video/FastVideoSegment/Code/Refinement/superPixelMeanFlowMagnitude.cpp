/** Function to compute the mean optical flow magnitude in each superpixel
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
			1 - Flow: a HxWx2 int16 matrix containing the optical flow
				between two pairs of frames 
			2 - Superpixel Map Previous: a HxW uint16 super pixel label map
		Optional inputs:
			3 - Number of superpixels in frame
*/

#include <matrix.h>
#include <mex.h>
#include <cmath>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[] )
{
	switch( nrhs )
	{
		case 2:
			if( mxGetClassID( prhs[ 0 ] ) != mxINT16_CLASS )
				mexErrMsgTxt( "Flow matrix must be of type int16" );
			if( mxGetClassID( prhs[ 1 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			break;
		case 3:
			if( mxGetClassID( prhs[ 0 ] ) != mxINT16_CLASS )
				mexErrMsgTxt( "Flow matrix must be of type int16" );
			if( mxGetClassID( prhs[ 1 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			if( mxGetClassID( prhs[ 2 ] ) != mxDOUBLE_CLASS )
				mexErrMsgTxt( "Number of superpixels must be of type double" );
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
	
	if( mxGetNumberOfDimensions( prhs[ 0 ] ) != 3 )
		mexErrMsgTxt( "Flow matrix must be 3-dimensional" );

	if( mxGetNumberOfDimensions( prhs[ 1 ] ) != 2 )
		mexErrMsgTxt( "Super pixel label map must be 2-dimensional" );

	int height = mxGetM( (mxArray *)prhs[ 1 ] );
	int width = mxGetN( (mxArray *)prhs[ 1 ] );
	//mexPrintf( "Matrix size: ( %i, %i )\n", height, width );
	
	short *flow = ( short * )mxGetData( prhs[ 0 ] );
	unsigned short *labelsMap = ( unsigned short * )mxGetData( prhs[ 1 ] );
	
	/** Get the number of superpixels in each frame */
	int labels;
	if( nrhs == 5 )
	{
		/** If number of superpixels is given, just use that */
		labels = ( int )( ( mxGetPr( prhs[ 2 ] )[ 0 ] ) );
	}
	else
	{
		/** If number of superpixels is not given, find the largest label */
		labels = 0;
		int point;
		for( int i = 0; i < height; i++ )
		{
			for( int j = 0; j < width; j++ )
			{
				point = j * height + i;
				if( labelsMap[ point ] > labels )
					labels = labelsMap[ point ];
			}
		}
	}
	
	/** Create the output (connectivity) matrix */
	mxArray *output =
		mxCreateNumericMatrix( labels, 1, mxSINGLE_CLASS, mxREAL );
	float *outData = ( float * )mxGetData( output );
	
	unsigned int *pixelCount = new unsigned int[ labels ];
	for( int i = 0; i < labels; i++ )
		pixelCount[ i ] = 0;
	
	/** Compute the connectivity matrix */
	int point;
	int label;
	int HxW = height * width;
	float magnitude;
	for( int i = 0; i < height; i++ )
	{
		for( int j = 0; j < width; j++ )
		{
			point = j * height + i;

			label = labelsMap[ point ] - 1;
			/** Make sure not to account for unlabeled pixels: label = 0 */
			if( ( label > -1 ) )
			{
				pixelCount[ label ]++;
				magnitude = sqrt( pow( (float)(flow[ point ]), 2.0 ) + pow( (float)(flow[ point + HxW ]), 2.0 ) );
				outData[ label ] += magnitude;
			}
		}
	}
	
	for( int i = 0; i < labels; i++ )
	{
		if( pixelCount[ i ] )
			outData[ i ] /= (float)( pixelCount[ i ] );
	}

	plhs[ 0 ] = output;

	delete [] pixelCount;

}
