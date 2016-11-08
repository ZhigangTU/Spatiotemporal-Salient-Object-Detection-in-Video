function [ MODEL, MAP ] = vectorQuantize( C, M, W )

% NEMA_VECTOR_QUANTIZE A vector quantization function that uses the
% binary split algorithm of Orchard and Bouman:
%
%  Color Quantization of Images, M. Orchard and C. Bouman, IEEE
%  Trans. on Signal Processing, Vol. 39, No. 12, pp, 2677--2690,
%  Dec. 1991.
%
%    [ MODEL, MAP ] = NEMA_COLOUR_QUANTIZE( C, M ) clusters the data
%    matrix C into M clusters.  Each column of the DxN matrix C is a
%    data point.  The output MODEL is a 1xN matrix that assigns
%    each data point to a cluster.  Each column i of MAP is the
%    mean of cluster i.
%
%    [ X, MAP ] = NEMA_COLOUR_QUANTIZE( I, M ) converts the RGB image
%    I into an indexed image X.  MAP contains at most M colours.
%
%    ... = NEMA_COLOUR_QUANTIZE( ..., W ) uses "subjectively weighted
%    TSE criteria" to quantize the data.

% Author: Nicholas Apostoloff <nema@robots.ox.ac.uk>
% Date: 23 May 05

    return_map = 0;
    if size( C, 3 ) > 1
        [ nrows, ncols, ndims ] = size( C );
        C = reshape( C, [], ndims )';
        C = im2double( C );
        return_map = 1;
    else
        [ndims, ncols] = size( C );
        C = double( C );
        nrows = 1;
    end
    npixels = nrows * ncols;

    if nargin < 3 || isempty( W )
      W = ones( 1, size( C, 2 ) );
    else
      W = W( : )';
    end % if nargin
    sumW = sum( W );

    if sumW
      W = W / sumW;
    else
      warning( 'Weights vector sums to zero.  Defaulting to ones' );
      W = ones( size( W ) ) / prod( size( W ) );
    end

    % remove data points with zero weight
    %idxw = find( W );
    %C = C( :, idxw );
    %W = W( idxw );

    if length( W ) ~= size( C, 2 )
      error( 'Weight vector does not match the data vector' );
    end

    % initial cluster
    [R0, m0] = ncq_Rm( C, W );
    N0 = sum( W );
    idx0 = ones( 1, size( C, 2 ) );
    [e0, lambda0] = ncq_e( C, R0, m0, W );

    for m = 2:M

      % find the node to split
      [t1, n] = max( lambda0 );
      en = e0( :, n );

      if ~t1 % if max( lambda0 ) is 0 then there is no variation left
        break;
      end

      % split node n
      nidx = find( idx0 == n );
      Cn = C( :, nidx );
      Wn = W( nidx );

      Nn = N0( n );
      mn = m0( :, n );
      qn = mn / Nn;
      Rn = R0( :, :, n );
      t1 = en' * Cn;
      nidx1 = find( t1 <= en' * qn );
      nidx2 = find( t1 > en' * qn );

      %  if ~length( nidx1 ) | ~length( nidx2 )
      %  if sum( Wn( nidx1 ) ) < 1 | sum( Wn( nidx2 ) ) < 1
      if length( nidx1 ) < 1 | length( nidx2 ) < 1
        break;
      end

      % update cluster parameters
      idx0( nidx( nidx2 ) ) = m;
      [ R0( :, :, n ), m0( :, n ) ] = ncq_Rm( Cn( :,  nidx1 ), Wn( nidx1 ) );
      [ e0( :, n ), lambda0( n ) ] = ncq_e( Cn( :, nidx1 ), R0( :, :, n ), ...
                        m0( :, n ), Wn( nidx1 ) );
      N0( n ) = sum( Wn( nidx1 ) );

      %[ R0( :, :, m ), m0( :, m ) ] = ncq_Rm( Cn( :,  nidx2 ), Wn( nidx2 ) );  
      %N0( m ) = sum( Wn( nidx2 ) );

      R0( :, :, m ) = Rn - R0( :, :, n );
      m0( :, m ) = mn - m0( :, n );
      N0( m ) = Nn - N0( n );

      [ e0( :, m ), lambda0( m ) ] = ncq_e( Cn( :, nidx2 ), R0( :, :, m ), ...
                        m0( :, m ), Wn( nidx2 ) );

    end % for m
    N0(N0==0) = 1;
    MAP = ( m0./repmat( N0, [ndims, 1] ) );
    if return_map
      MAP = MAP';
    end
    MODEL = reshape( idx0, nrows, ncols );

end

function [R, m] = ncq_Rm( C, W )
C2 = repmat( W, size( C, 1 ), 1 ) .* C;
R = C2 * C';
m = sum( C2, 2 );

end

function [en, lambda] = ncq_e( Cn, Rn, mn, Wn )  
Nn = sum( Wn );

qn = mn / Nn;

Rhat = Rn - mn * mn' / Nn;

% Rhat should be real symmetric - make is so
%  this should stop a bug where complex eigenvalues where returned
Rhat = (Rhat + Rhat') * 0.5;

[ V, D ] = eigs( Rhat );
if any( ~isreal( V( :, 1 ) ) )
  warning( 'Complex eigenvectors found in Rhat' );
end
en = V( :, 1 );

%%%% In their paper they do not specify that the variance of each
% mode should be calculated using the weighted points, however I
% think that this is incorrect.

% Their method
%lambda = sum((en' * (Cn - repmat( qn, [1, size( Cn, 2 )] ))).^2);

% My method
%lambda = sum((en' * (Cn - repmat( qn, [1, size( Cn, 2 )] ))).^2 .* Wn)

% However, the variance of the data points along an eigenvector is
% equal to the eigenvalue of that eigenvector
lambda = D( 1 );

end