% Function to auto-calibrate the motion bounds lambda weight
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

function bestEdges = calibrateProbabilityEdge( flowframe, lambda )

    epsilon = 0.1;

    % Starting edge sensitivity point
    if( ~exist( 'lambda', 'var' ) || isempty( lambda ) )
        lambda = 0.7;
    end
    
    [ height, width, ~ ] = size( flowframe );
    
    magnitude = getMagnitude( getFlowGradient( flowframe ) );
    rotBoundary = getFlowDifference( flowframe );
    
    bestQuality = -1;
    bestEdges = zeros( height, width );
    previousEdges = false( height, width );
    for i = lambda: epsilon: 1.5 

        edges = getEdge( magnitude, rotBoundary, i );
        thresholdedEdges = edges > 0.5;
        
        if( thresholdedEdges == previousEdges )
            continue;
        else
            inPoints = getInPoints( thresholdedEdges, [], true );
            quality = getFrameQuality( inPoints > 4 );
            previousEdges = thresholdedEdges;
        end
        
        if( quality > bestQuality )
            bestQuality = quality;
            bestEdges = edges;
        end
    end
   
end

function result = getEdge( magnitude, rotBoundary, gradLambda )
    gradBoundary = 1 - exp( -gradLambda * magnitude );        

    large = gradBoundary > 0.6 ;
    medium = gradBoundary <= 0.6 & gradBoundary > 0.25;
    result = 0.1 * gradBoundary;
    result(large) = gradBoundary(large);
    result(medium) = (gradBoundary(medium) .* rotBoundary(medium));
end
