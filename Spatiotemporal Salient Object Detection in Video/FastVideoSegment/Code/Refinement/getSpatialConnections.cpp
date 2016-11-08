/** Function to compute the spatial connection in each frame
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
			1 - Superpixels: a cell array of length frames. Each cell i contains a
				HxW uint32 superpixel superpixel map of frame i.
			2 - Number of superpixels: a double value containing the total number of
				superpixel superpixels.

		Outputs:
			1 - Source superpixels:
			2 - Target superpixels:
*/

#include <matrix.h>
#include <mex.h>

#include <set>

//#define DEBUG_MODE

#define USAGE_NOTE " "

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	unsigned int frames, height, width, point, pointNext, superpixel, superpixelNext, count;
    unsigned long long superpixels, connection;
	unsigned int *superpixelMap, *sources, *targets;

	std::set<unsigned long long> connections;

	if( nrhs == 2 )
	{
		/** Assert that inputs are cell arrays */
		if( !mxIsCell( prhs[ 0 ] ) )
			mexErrMsgTxt( USAGE_NOTE );

		/** Assert that input cell arrays are not empty */
		frames = mxGetNumberOfElements( prhs[ 0 ] );
		if( frames == 0 )
			mexErrMsgTxt( USAGE_NOTE );

		/** Assert that cell contents are of correct data type */
		for( int frame = 0; frame < frames; frame++ )
		{
			if( mxGetClassID( mxGetCell( prhs[ 0 ], frame ) ) != mxUINT32_CLASS )
				mexErrMsgTxt( USAGE_NOTE );
		}
		if( mxGetClassID( prhs[ 1 ] ) != mxDOUBLE_CLASS )
				mexErrMsgTxt( USAGE_NOTE );
	}
	else
	{
			mexErrMsgTxt( USAGE_NOTE );
	}

	switch( nlhs )
	{
		case 0:
			break;
		case 2:
			break;
		default:
			mexErrMsgTxt( USAGE_NOTE );
	}

	height = mxGetM( mxGetCell( prhs[ 0 ], 0 ) );
	width = mxGetN( mxGetCell( prhs[ 0 ], 0 ) );
	superpixels = ( unsigned long long )( ( double * )mxGetData( prhs[ 1 ] ) )[ 0 ];
	#ifdef DEBUG_MODE
	mexPrintf( "getSpatialConnections: Number of frames: %i\n", frames );
	mexPrintf( "getSpatialConnections: Frame size: %ix%i\n", height, width );
	mexPrintf( "getSpatialConnections: Total number of superpixels: %i\n", superpixels );
	mexEvalString( "pause(0.001)" );
	#endif

	for( unsigned int frame = 0; frame < frames; frame++ )
	{
		superpixelMap = ( unsigned int * )( mxGetData( mxGetCell( prhs[ 0 ], frame ) ) );

		for( unsigned int i = 0; i < height; i++ )
		{
			for( unsigned int j = 0; j < width; j++ )
			{
				point = j * height + i;
				superpixel = superpixelMap[ point ] - 1;

				if( superpixel < 0 && superpixel >= superpixels )
					mexErrMsgTxt( "getSpatialConnections: Superpixels found outside given range" );

				if( i > 0 )
				{
					pointNext = j * height + ( i - 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( j > 0 )
				{
					pointNext = ( j - 1 ) * height + i;
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( i < height - 1 )
				{
					pointNext = j * height + ( i + 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( j < width - 1 )
				{
					pointNext = ( j + 1 ) * height + i;
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( i > 0 && j > 0 )
				{
					pointNext = ( j - 1 ) * height + ( i - 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( i > 0 && j < width - 1 )
				{
					pointNext = ( j + 1 ) * height + ( i - 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( i < height - 1 && j > 0 )
				{
					pointNext = ( j - 1 ) * height + ( i + 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}

				if( i < height - 1 && j < width - 1 )
				{
					pointNext = ( j + 1 ) * height + ( i + 1 );
					superpixelNext = superpixelMap[ pointNext ] - 1;

					if( superpixel < superpixelNext )
						connections.insert( superpixel + superpixels * superpixelNext );
				}
			}
		}
	}

	#ifdef DEBUG_MODE
	mexPrintf( "getSpatialConnections: Number of connections: %i\n", connections.size() );
	mexEvalString( "pause(0.001)" );
	#endif

	#ifdef DEBUG_MODE
	mexPrintf( "getSpatialConnections: Allocating memory for outputs...\n" );
	mexEvalString( "pause(0.001)" );
	#endif

	mxArray *sourcesMxArray = mxCreateNumericMatrix( connections.size(), 1, mxUINT32_CLASS, mxREAL );
	mxArray *targetsMxArray = mxCreateNumericMatrix( connections.size(), 1, mxUINT32_CLASS, mxREAL );

	sources = ( unsigned int * )mxGetData( sourcesMxArray );
	targets = ( unsigned int * )mxGetData( targetsMxArray );

	count = 0;
	for( std::set<unsigned long long>::iterator i = connections.begin(); i != connections.end(); i++ )
	{
		connection = *i;

		targets[ count ] = connection / superpixels;
		sources[ count ] = connection % superpixels;

		count++;
	}

	plhs[ 0 ] = sourcesMxArray;
	plhs[ 1 ] = targetsMxArray;

}

