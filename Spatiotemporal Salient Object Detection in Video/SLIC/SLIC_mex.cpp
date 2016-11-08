/** A wrapper function for the SLIC superpixel oversegmentation method
% 
% Copyright (c) 2013, Anestis Papazoglou
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met: 
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer. 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


#include <string>
#include <vector>
#include <algorithm>

#include <matrix.h>
#include <mex.h>

#include "SLIC.h"

#define OPAQUE 0xFF0000
#define SQUARE255 65025

#define USAGE_NOTE "USAGE: \n"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	bool isRGB;
	unsigned char R, G, B;
	unsigned int height, width, pixels, pixelsX2, dimensions;
	int superpixelsNumber, labelsNumber;
	double compactness;
	unsigned char *image;
	unsigned short *output;
	unsigned int *aRGBimage;
	int *superpixels;

	const mxArray *imageMxArray;
	const mxArray *superpixelsNumberMxArray;
	mxArray *superpixelsMxArray;

	if( nrhs == 2 || nrhs == 3 )
	{
		imageMxArray = prhs[ 0 ];
		superpixelsNumberMxArray = prhs[ 1 ];

		if( mxGetClassID( imageMxArray ) != mxUINT8_CLASS )
				mexErrMsgTxt( USAGE_NOTE );
				
		if( mxGetClassID( superpixelsNumberMxArray ) != mxDOUBLE_CLASS )
				mexErrMsgTxt( USAGE_NOTE );

		dimensions = mxGetNumberOfDimensions( imageMxArray );

		if( dimensions == 3 )
		{
			//TODO: make sure 3rd dimensions contains 3 colours
			isRGB = true;
			height = mxGetM( imageMxArray );
			width = mxGetN( imageMxArray ) / 3;
		}
		else if( dimensions == 1 )
		{
			isRGB = false;
			height = mxGetM( imageMxArray );
			width = mxGetN( imageMxArray );
		}
		else
		{
			mexErrMsgTxt( USAGE_NOTE );
		}
		pixels = height * width;
		pixelsX2 = 2 * pixels;
		
		image = ( unsigned char * )mxGetData( imageMxArray );
		superpixelsNumber = ( unsigned int )( ( double * )mxGetData( superpixelsNumberMxArray ) )[ 0 ];

		if( ( superpixelsNumber < 1 ) || ( superpixelsNumber >= pixels ) )
			mexErrMsgTxt( USAGE_NOTE );
		
		if( nrhs == 3 )
		{
			if( mxGetClassID( prhs[ 2 ] ) == mxDOUBLE_CLASS	 )
				compactness = ( ( double * )mxGetData( prhs[ 2 ] ) )[ 0 ];
			else
				mexErrMsgTxt( USAGE_NOTE );
		}
		else
		{
			compactness = 20;
		}
	}
	else
	{
		mexErrMsgTxt( USAGE_NOTE );
	}
	
	superpixelsMxArray = mxCreateNumericMatrix( height, width, mxUINT16_CLASS, mxREAL ); 
	//superpixels = ( int * )mxGetData( superpixelsMxArray );
	aRGBimage = new unsigned int[ pixels ];
	
	for( unsigned int i = 0; i < pixels; i++ )
	{
		R = image[ i ];
		G = image[ i + pixels ];
		B = image[ i + pixelsX2 ];
		
		aRGBimage[ i ] = OPAQUE + SQUARE255 * R + 255 * G + B;
	}

	SLIC slic;

	slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels( aRGBimage, height, width, superpixels, labelsNumber, superpixelsNumber, compactness );
	//mexEvalString( "pause" );
	//mxSetData( superpixelsMxArray, ( void * )superpixels );
	//mxSetM( superpixelsMxArray, height );
	//mxSetN( superpixelsMxArray, width );
	
	output = ( unsigned short * )mxGetData( superpixelsMxArray );
	for( int i = 0; i < pixels; i++ )
	{
		output[ i ] = ( unsigned short )( superpixels[ i ] + 1 );
	}
	
	plhs[ 0 ] = superpixelsMxArray;
	
	delete [] aRGBimage;
	
	return;
}

