% % Function to load all frames in a shot
% %
% %    Copyright (C) 2013  Anestis Papazoglou
% %
% %    You can redistribute and/or modify this software for non-commercial use
% %    under the terms of the GNU General Public License as published by
% %    the Free Software Foundation, either version 3 of the License, or
% %    (at your option) any later version.
% %
% %    This program is distributed in the hope that it will be useful,
% %    but WITHOUT ANY WARRANTY; without even the implied warranty of
% %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% %    GNU General Public License for more details.
% %
% %    You should have received a copy of the GNU General Public License
% %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% %
% %    For commercial use, contact the author for licensing options.
% %
% %    Contact: a.papazoglou@sms.ed.ac.uk
% 
% function frames = readAllFrames( options, range )
% 
%     start = options.ranges( range );
%     stop = options.ranges( range + 1 ) - 1;
% 
%     frames = cell( stop - start + 1, 1 );
%     for i = start: stop
%         index = i - start + 1;
%         frames{ index } = readFrame( options, i );
%     end
% 
% end


function  [frames, names, height, width, nframe ] = readAllFrames(folderName)
    
    frameFiles = imdir(fullfile(folderName));
    nframe = length(frameFiles);
    frames = cell(nframe, 1);
    names = cell(nframe, 1);
    for index = 1: nframe 
        [~, frameName] = fileparts(frameFiles(index).name);
        if exist(fullfile(folderName, [frameName '.png']),'file')
            frame = uint8(imread(fullfile(folderName, [frameName '.png'])));
        elseif exist(fullfile(folderName, [frameName '.jpg']),'file')
            frame = uint8(imread(fullfile(folderName, [frameName '.jpg'])));
        elseif exist(fullfile(folderName, [frameName '.bmp']),'file')
            frame = uint8(imread(fullfile(folderName, [frameName '.bmp'])));
        end
        
%         frame=imfilter(frame,fspecial('gaussian',7,1.5),'same','replicate');

        frames{ index } = uint8(frame);
        names{ index } =  frameName;
    end
    [ height,width ] = size(frame(:,:,1));
end
