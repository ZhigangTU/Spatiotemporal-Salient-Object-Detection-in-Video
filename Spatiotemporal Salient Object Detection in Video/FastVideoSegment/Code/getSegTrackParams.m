% Function to provide parameters for SegTrack videos. Valid video names are
%   'birdfall', 'cheetah', 'girl', 'monkey' and 'parachute'
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

function params = getSegTrackParams( videoid )

    params.fadeout = 0.0001;
    params.foregroundMixtures = 5;
    params.backgroundMixtures = 8;
    params.maxIterations = 1;
    
    params.locationNorm = 0.95;
    
    if(strcmp(videoid, 'birdfall2') || strcmp(videoid, 'girl')) %|| strcmp( videoid, 'monkey' )
    
        params.locationWeight = 1.5;
        params.spatialWeight  = 250;
        params.temporalWeight = 1000;
    elseif(strcmp(videoid, 'cheetah') || strcmp(videoid, 'parachute'))       
        params.locationWeight = 1;
        params.spatialWeight  = 2000;
        params.temporalWeight = 1000;
    elseif(strcmp(videoid, 'monkeydog'))
        params.locationWeight = 5;
        params.temporalWeight = 4000; % Test param
        params.spatialWeight  = 5000; % Test param
    else
        params = getDefaultParams();
%         error( 'getSegTrackParams: "%s" is not a valid videoid', videoid );
    end

end
