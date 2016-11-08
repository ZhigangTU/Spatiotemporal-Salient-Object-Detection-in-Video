/** Function to compute the superpixel connections between two subsequent frames
		based on the optical flow
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
				for the first frame
			3 - Superpixel Map Next: a HxW uint16 super pixel label map
				for the second frame
		Optional inputs:
			4 - Number of superpixels in first frame
			5 - Number of superpixels in second frame
			
*/

#include <matrix.h>
#include <mex.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[] )
{
	switch( nrhs )
	{
		case 3:
			if( mxGetClassID( prhs[ 0 ] ) != mxINT16_CLASS )
				mexErrMsgTxt( "Flow matrix must be of type int16" );
			if( mxGetClassID( prhs[ 1 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			if( mxGetClassID( prhs[ 2 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			break;
		case 5:
			if( mxGetClassID( prhs[ 0 ] ) != mxINT16_CLASS )
				mexErrMsgTxt( "Flow matrix must be of type int16" );
			if( mxGetClassID( prhs[ 1 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			if( mxGetClassID( prhs[ 2 ] ) != mxUINT16_CLASS )
				mexErrMsgTxt( "Super pixel label map must be of type uint16" );
			if( mxGetClassID( prhs[ 3 ] ) != mxDOUBLE_CLASS )
				mexErrMsgTxt( "Number of superpixels must be of type double" );
			if( mxGetClassID( prhs[ 4 ] ) != mxDOUBLE_CLASS )
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
		
	if( mxGetNumberOfDimensions( prhs[ 2 ] ) != 2 )
		mexErrMsgTxt( "Super pixel label map must be 2-dimensional" );

	int height = mxGetM( (mxArray *)prhs[ 1 ] );
	int width = mxGetN( (mxArray *)prhs[ 1 ] );
	//mexPrintf( "Matrix size: ( %i, %i )\n", height, width );
	
	short *flow = ( short * )mxGetData( prhs[ 0 ] );
	unsigned short *previousFrame = ( unsigned short * )mxGetData( prhs[ 1 ] );
	unsigned short *nextFrame = ( unsigned short * )mxGetData( prhs[ 2 ] );
	
	/** Get the number of superpixels in each frame */
	int labelsPrevious;
	int labelsNext;
	if( nrhs == 5 )
	{
		/** If number of superpixels is given, just use that */
		labelsPrevious = ( int )( ( mxGetPr( prhs[ 3 ] )[ 0 ] ) );
		labelsNext = ( int )( ( mxGetPr( prhs[ 4 ] )[ 0 ] ) );
	}
	else
	{
		/** If number of superpixels is not given, find the largest label */
		labelsPrevious = 0;
		labelsNext = 0;
		int point;
		for( int i = 0; i < height; i++ )
		{
			for( int j = 0; j < width; j++ )
			{
				point = j * height + i;
				if( previousFrame[ point ] > labelsPrevious )
					labelsPrevious = previousFrame[ point ];
				
				if( nextFrame[ point ] > labelsNext )
					labelsNext = nextFrame[ point ];
			}
		}
	}
	
	/** Create the output (connectivity) matrix */
	mxArray *output =
		mxCreateNumericMatrix( labelsPrevious, labelsNext, mxUINT32_CLASS, mxREAL );
	unsigned int *outData = ( unsigned int * )mxGetData( output );
	
	/** Compute the connectivity matrix */
	int point, nextPoint;
	int xNext, yNext;
	int previousLabel, nextLabel;
	int HxW = height * width;
	for( int i = 0; i < height; i++ )
	{
		for( int j = 0; j < width; j++ )
		{
			point = j * height + i;
			xNext = i + flow[ point ];
			yNext = j + flow[ point + HxW ];
			if( xNext > -1 && xNext < height &&
				yNext > -1 && yNext < width )
			{
				nextPoint = yNext * height + xNext;
				previousLabel = previousFrame[ point ] - 1;
				nextLabel = nextFrame[ nextPoint ] - 1;
				/** Make sure not to account for unlabeled pixels: label = 0 */
				if( ( previousLabel > -1 ) && ( nextLabel > -1 ) )
					outData[ previousLabel + nextLabel * labelsPrevious ]++;
			} 
		}
	}

	plhs[ 0 ] = output;

}
