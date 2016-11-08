/** Function to add the weights of each unique colour code
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

#include <matrix.h>
#include <mex.h>

#include <map>

#define USAGE_NOTIFICATION "Error: incorrect function call\n\n" \
	"USAGE:\t[ colours, weights ] = findUniqueColourWeights( colours, weights )\n\n" \
	"\tcolours: a Nx3 uint8 matrix containing the colour list\n" \
	"\tweights: a Nx1 single real valued array containing the corresponding weight " \
	"of each colour\n"

#define RANGE 256
#define DOUBLE_RANGE 65536

void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[] )
{
	if( nrhs != 2 )
		mexErrMsgTxt( USAGE_NOTIFICATION );
		
	if( nlhs > 2 )
		mexErrMsgTxt( USAGE_NOTIFICATION );

	if( mxGetClassID( prhs[ 0 ] ) != mxUINT8_CLASS )
		mexErrMsgTxt( USAGE_NOTIFICATION );
	if( mxGetClassID( prhs[ 1 ] ) != mxSINGLE_CLASS )
		mexErrMsgTxt( USAGE_NOTIFICATION );

	unsigned int samples = mxGetM( prhs[ 0 ] );
	unsigned int dimensions = mxGetN( prhs[ 0 ] );
	unsigned int doubleSamples = 2 * samples;
	
	if( samples != mxGetM( prhs[ 1 ] ) ||	dimensions!= 3 || mxGetN( prhs[ 1 ] ) != 1 )
		mexErrMsgTxt( USAGE_NOTIFICATION );
	
	unsigned char *colours = ( unsigned char * )mxGetData( prhs[ 0 ] );
	float *weights = ( float * )mxGetData( prhs[ 1 ] );
	
	//typedef std::tr1:unordered_map< int, double > HashMap;
	typedef std::map< unsigned int, float > ColourMap;
	ColourMap colourMap;

	unsigned int key;
	for( int i = 0; i < samples; i++ )
	{
		key = colours[ i ] + RANGE * colours[ i + samples ] + DOUBLE_RANGE * colours[ i + doubleSamples ];
		colourMap[ key ] += weights[ i ];
	}

	unsigned int unique = colourMap.size();
	unsigned int doubleUnique = 2 * unique;
	
	mxArray *uniqueColoursMxArray =	mxCreateNumericMatrix( unique, 3, mxUINT8_CLASS, mxREAL );
	unsigned char *uniqueColours = ( unsigned char * )mxGetData( uniqueColoursMxArray );
	
	mxArray *weightSumsMxArray = mxCreateNumericMatrix( unique, 1, mxSINGLE_CLASS, mxREAL );
	float *weightSums = ( float * )mxGetData( weightSumsMxArray );
	
	ColourMap::const_iterator iterator;
	unsigned int count = 0;
	unsigned char red, green, blue;
	for( iterator = colourMap.begin(); iterator != colourMap.end(); iterator++ )
	{
		key = ( *iterator ).first;
		blue = key / DOUBLE_RANGE;
		green = ( key - DOUBLE_RANGE * blue ) / RANGE;
		red = key - RANGE * green - DOUBLE_RANGE * blue;
		
		uniqueColours[ count ] = red;
		uniqueColours[ count + unique ] = green;
		uniqueColours[ count + doubleUnique ] = blue;

		weightSums[ count ] = ( *iterator ).second;

		count++;
	}

	plhs[ 0 ] = uniqueColoursMxArray;
	plhs[ 1 ] = weightSumsMxArray;

}
