% Script to compile all required sofware
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
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Compile Matlab wrappers for external code
cd Code/External/maxflow/
mex maxflow_mex_optimisedWrapper.cpp maxflow-v3.0/graph.cpp maxflow-v3.0/maxflow.cpp

cd ../SLIC/
mex SLIC_mex.cpp SLIC.cpp

cd ../sundaramECCV2010/
mex sundaramECCV10_ldof_GPU_mex.cpp NMath.cpp -lldof_gpu

% Build TurboPixels
cd ../TurboPixels
make

% Compile code for VideoRapidSegment
cd ../../Integral' I'ntersections
mex integralIntersections.cpp

cd ../Motion' B'ounds
mex cleanBoundaries.cpp

cd ../Refinement
mex findUniqueColourWeights.cpp
mex getSpatialConnections.cpp
mex getSuperpixelStats.cpp
mex getTemporalConnections.cpp
mex superPixelConnectivity.cpp
mex superpixelInRatio.cpp
mex superPixelMeanFlowMagnitude.cpp

cd ../Utilities
mex fillTurboPixels.cpp

cd ../..
