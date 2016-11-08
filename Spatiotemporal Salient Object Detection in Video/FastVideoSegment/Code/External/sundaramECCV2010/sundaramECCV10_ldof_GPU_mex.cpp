/** Wrapper function for the sundaramECCV10 GPU optical flow estimation method
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
#include <ctime>
#include "CTensorMatlab.h"
#include "CFilter.h"
#include "ldof.h"

#include <matrix.h>
#include <mex.h>

#define USAGE_NOTE "USAGE: \n"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	unsigned int height, width, dimensions;

	const mxArray *prevFrameMxArray;
	const mxArray *nextFrameMxArray;
	
	mxArray *forwardFlowMxArray;
	mxArray *backwardFlowMxArray;
	
	if( nrhs == 2 )
	{
		prevFrameMxArray = prhs[ 0 ];
		nextFrameMxArray = prhs[ 1 ];
		
		if( mxGetClassID( prevFrameMxArray ) != mxUINT8_CLASS ||
			mxGetClassID( nextFrameMxArray ) != mxUINT8_CLASS )
		{
			mexErrMsgTxt( USAGE_NOTE );
		}
		
		dimensions = mxGetNumberOfDimensions( prevFrameMxArray );
		
		if( dimensions != mxGetNumberOfDimensions( nextFrameMxArray ) )
			mexErrMsgTxt( USAGE_NOTE );

		if( dimensions != 3 )
			mexErrMsgTxt( USAGE_NOTE );
		
		height = mxGetM( prevFrameMxArray );
		width = mxGetN( prevFrameMxArray ) / 3;
		
		if( height != mxGetM( prevFrameMxArray ) || width != mxGetN( prevFrameMxArray ) / 3 )
			mexErrMsgTxt( USAGE_NOTE );
	}
	else
	{
		mexErrMsgTxt( USAGE_NOTE );
	}

	CTensorMatlab<float>* prevFrame = new CTensorMatlab<float>;
	CTensorMatlab<float>* nextFrame = new CTensorMatlab<float>;
	
	prevFrame->copyFromMxArrayImage( prevFrameMxArray );
	nextFrame->copyFromMxArrayImage( nextFrameMxArray );
	
	NFilter::recursiveSmoothX( *prevFrame, 0.8f );
	NFilter::recursiveSmoothY( *prevFrame, 0.8f );
	
	NFilter::recursiveSmoothX( *nextFrame, 0.8f );
	NFilter::recursiveSmoothY( *nextFrame, 0.8f );
	
	CTensorMatlab<float> fflow;
	CTensorMatlab<float> bflow;
	
	ldof( *prevFrame, *nextFrame, fflow, bflow );
	
	const mwSize flowSize[] = { height, width, 2 }; 
	
	forwardFlowMxArray = mxCreateNumericArray( 3, flowSize, mxSINGLE_CLASS, mxREAL );
	backwardFlowMxArray = mxCreateNumericArray( 3, flowSize, mxSINGLE_CLASS, mxREAL );
	
	fflow.copyToMxArray( forwardFlowMxArray );
	bflow.copyToMxArray( backwardFlowMxArray );
	
	plhs[ 0 ] = forwardFlowMxArray;
	plhs[ 1 ] = backwardFlowMxArray;
	
	delete prevFrame;
	delete nextFrame;
}
