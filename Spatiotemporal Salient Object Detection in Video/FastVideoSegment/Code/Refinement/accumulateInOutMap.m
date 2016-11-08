% Function to accumulate the inside rations based on the optical flow
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
%    Contact: a.papazoglou@sms.ed.ac.uk

function [ result, tags ] = accumulateInOutMap( params, data )
    
    if( exist( 'params', 'var' ) && isfield( params, 'lambda' ) )
        lambda = params.lambda;
    else
        lambda = 2; %5
    end
    
    flow = data.flow;
    superpixels = data.superpixels;
    map = data.inMaps;
    connectivity = getSuperpixelConnectivity( flow, superpixels );

    spixelRatio = getSuperpixelInRatio( superpixels, map );
    
    frames = length( flow );
    
    % fprintf( 'Propagating metric forward...\n' );
    output = cell( frames, 1 );
    foreGain = cell( frames, 1 );
    output{ 1 } = double( spixelRatio{ 1 }' );
    foreGain{ 1 } = zeros( length( spixelRatio{ 1 } ), 1 );
    for( frame = 1: frames - 1 )

        % Find number of connections between superpixels
        connect = connectivity{ frame };
        
        % Find the ratio of connections leading to a specific superpixel
        normFactor =  1 ./ sum( connect, 2 );
        % Check for inf - Turbopixels may return empty superpixel sets
        normFactor( isinf( normFactor ) ) = 0;
        normFactor = sparse( 1: length( normFactor ), ...
            1: length( normFactor ), normFactor );
        connect = normFactor * connect;
        
        alpha = superPixelMeanFlowMagnitude( int16( 100 * ...
            getFlowGradient( flow{ frame } ) ), ...
            superpixels{ frame } ) / 100;
        % This should be normalised for frame size
        alpha = double( exp( - lambda * alpha ) );
        alpha = sparse( 1: length( alpha ), 1: length( alpha ), alpha );

        alphaConnect = alpha * connect;
        sumAlphaConnect = sum( alphaConnect );
        foreGain{ frame + 1 } = ( output{ frame } * alphaConnect ) ./ sumAlphaConnect;
        foreGain{ frame + 1 }( sumAlphaConnect == 0 ) = 0;
        foreGain{ frame + 1 } = foreGain{ frame + 1 }';
        output{ frame + 1 } = ( foreGain{ frame + 1 } + ...
            double( spixelRatio{ frame + 1 } ) )';

    end

    % fprintf( 'Propagating metric backward...\n' );
    backGain = cell( frames, 1 );
    output = cell( frames, 1 );
    output{ frames } = double( spixelRatio{ frames } );
    backGain{ frames } = zeros( length( spixelRatio{ frames } ), 1 );
    for( frame = frames - 1: -1: 1 )

        connect = connectivity{ frame };

        normFactor = 1 ./ sum( connect, 1 );
        % Check for inf - Turbopixels may return empty superpixel sets
        normFactor( isinf( normFactor ) ) = 0;
        normFactor = sparse( 1: length( normFactor ), ...
            1: length( normFactor ), normFactor );
        connect = connect * normFactor;

        alpha = superPixelMeanFlowMagnitude( int16( 100 * ...
            getFlowGradient( flow{ frame } ) ), ...
            superpixels{ frame } ) / 100;
        alpha = double( exp( - lambda * alpha ) );
        alpha = sparse( 1: length( alpha ), 1: length( alpha ), alpha );

        alphaConnect = alpha * connect;
        sumAlphaConnect = sum( alphaConnect, 2 );
        backGain{ frame } = ( alphaConnect * output{ frame + 1 } ) ./ sumAlphaConnect;
        backGain{ frame }( sumAlphaConnect == 0 ) = 0;

        output{ frame } = backGain{ frame } + double( spixelRatio{ frame } );

    end

    result = cell( frames, 1 );
    tags = cell( frames, 1 );
    for( i = 1: frames )
        tags{ i } = spixelRatio{ i } + foreGain{ i } + backGain{ i };
        
        result{ i } = tags{ i }( superpixels{ i } );
    end
    
end
