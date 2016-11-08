function [chnsReg,chnsSim] = edgesChns( I, opts )
% Compute features for structured edge detection.
%
% For an introductory tutorial please see edgesDemo.m.
%
% USAGE
%  [chnsReg,chnsSim] = edgesChns( I, opts )
%
% INPUTS
%  I          - [h x w x nchns] input data
%  opts       - structured edge model options
%
% OUTPUTS
%  chnsReg    - [h x w x nChannel] regular output channels
%  chnsSim    - [h x w x nChannel] self-similarity output channels
%
% Code written by Philippe Weinzaepfel, 2015, 
% for motion boundaries detection
% Licensed under the MSR-LA Full Rights License [see license.txt]

shrink=opts.shrink; chns=cell(1,opts.nChns); k=0;

mode = opts.mode;

mode_nchannels = [3,2,2,2,2]; % nchannels for image, flow, warping, bwd-flow, bwd-warping
mode_gradient = [2,1,0,1,0]; % gradient mode (0 => no, 1 => coarse scale, 2 => coarse and fine scales)
if strcmp(mode, 'Color'),
    assert( size(I,3)==3 );
    nmodes = 1;
elseif strcmp(mode, 'Color+Flow');
    assert( size(I,3)==5 );
    nmodes = 2;
elseif strcmp(mode, 'Color+Flow+Warping');
    assert( size(I,3)==7 );
    nmodes = 3;
elseif strcmp(mode, 'Color+Flow+Warping+Backward');
    assert( size(I,3)==11 );
    nmodes = 5;
else
    assert( 1==0 );
end

Is = I;
currentchannel=1;

for t=1:nmodes,
  I = Is(:,:,currentchannel:(currentchannel+mode_nchannels(t)-1)); 
  currentchannel = currentchannel+mode_nchannels(t);
  if(size(I,3)==3), I=rgbConvert(I,'luv'); end % 3 => color image
  Ishrink=imResample(I,1/shrink); k=k+1; chns{k}=Ishrink;
  if(mode_gradient(t)>0), % 0 no gradient
      for i = (3-mode_gradient(t)):2, % 1 coarse scale only, 2 coarse+fine
        s=2^(i-1);
        if(s==shrink), I1=Ishrink; else I1=imResample(I,1/s); end
        I1 = convTri( I1, opts.grdSmooth );
        if( size(I,3)==2 ), % flow
            [Mx,Ox] = gradientMag( I1, 1, opts.normRad, .01 );
            [My,Oy] = gradientMag( I1, 2, opts.normRad, .01 );
            M = sqrt(Mx.^2+My.^2);
            Ox = Ox*2;
            Oy = Oy*2;
            Mx = Mx+1e-6;
            My = My+1e-6;
            cosO = (Mx.*cos(Ox) + My.*cos(Oy)) ./ (Mx+My);
            sinO = (Mx.*sin(Ox) + My.*sin(Oy)) ./ (Mx+My);
            O = (atan2(sinO, cosO)+pi)/2;
        else % image
            [M,O] = gradientMag( I1, 0, opts.normRad, .01 );
        end   
        H = gradientHist( M, O, max(1,shrink/s), opts.nOrients, 0 );
        k=k+1; chns{k}=imResample(M,s/shrink);
        k=k+1; chns{k}=imResample(H,max(1,s/shrink));
      end
  end
end
chns=cat(3,chns{1:k}); assert(size(chns,3)==opts.nChns);
chnSm=opts.chnSmooth/shrink; if(chnSm>1), chnSm=round(chnSm); end
simSm=opts.simSmooth/shrink; if(simSm>1), simSm=round(simSm); end
chnsReg=convTri(chns,chnSm); chnsSim=convTri(chns,simSm);

end
