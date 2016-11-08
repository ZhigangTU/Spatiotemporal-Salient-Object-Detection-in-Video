function [ux,uy,Gm] = findGaCenter(Pb_map, SegCoord)
 
% load ObjectsPbWeight.mat
% load TargetObjectsBoxs.mat
% load TargetObjectsLabels.mat
% load SpixelList.mat
% noFrameImg = imread('778-62_l146.jpg');

% Debugging Code
%% SuperPixel Location
% Pb_Temp = frameoptFgConSacyMaps{1,1};
% SegCoord = cell2mat(frameTargetBox);
% pixelList = framePixelList{1,1};
% [h,w,d] = size(noFrameImg);
% idxImg = zeros(h,w);
% for idx = 1 : size(pixelList,1)
%     temp = cell2mat(pixelList(idx));
%     for spNum = 1 : size(temp,1)
%         idxImg(temp(spNum)) = idx;
%     end
% end

%% Number of BBoxes within a frame
idx_b = find(SegCoord(:,1) == SegCoord(1,1));
numOfBBox = size(idx_b,1);

%% Super Pixel Distribution within a BBox
Pb_Temp = Pb_map;

% pixelList = Spixelist;
% [h,w,d] = size(OriImg);
% idxImg = zeros(h,w);
% for idx = 1 : size(pixelList,1)
%     temp = cell2mat(pixelList(idx));
%     for spNum = 1 : size(temp,1)
%         idxImg(temp(spNum)) = idx;
%     end
% end


% idxSp = cellfun(@(x) [rem(x,720),round(x/720)+1],pixelList,'UniformOutput',false);
% idxSp = cellfun(@(x) [x(idx,:),x],idxSp,'UniformOutput',false);
% SpPb = load('SpixelPbWeight.mat');
% SpPb_temp = SpPb.frameoptwCtr{1, 1};

ux = cell(numOfBBox,1);
uy = cell(numOfBBox,1);
Gm = cell(numOfBBox,1);

%% For each Bounding Box
for i = 1 : numOfBBox;
    % Tailor the BBox data 
    TempCoord = SegCoord(idx_b(i),:);
    TempCoord = TempCoord(2:5);
    % Transfer to x range and y range
    % RegularCoord = [TempCoord(2),TempCoord(1),TempCoord(4)-TempCoord(2),TempCoord(3)-TempCoord(1)];
    XCoord = TempCoord(:,2) : TempCoord(:,4);
    YCoord = TempCoord(:,1) : TempCoord(:,3);
    
    %Spread the X,Y Coordinates
    F_XCoord = repmat(XCoord,size(YCoord,2),1);
    F_YCoord = repmat(YCoord',1,size(XCoord,2));
    
    % Pb within the region of interest
    Pb = Pb_Temp(YCoord,XCoord);

    
%     Spidx = idxImg(YCoord,XCoord);
%     Sp_unique = unique(Spidx);
%     num_Sp = histc(reshape(Spidx,size(Spidx,1)*size(Spidx,2),1),Sp_unique);
%     ref_Sp = SpPb_temp(Sp_unique);
%     Sp_pb = [Sp_unique,ref_Sp,num_Sp];

    % Load the Saliency map
    % SegLable = cell2mat(frameTargetLabels(end,1));
    % SegLable = cell2mat(frameTargetLabels);
    %imshow(SegLable)

    eta = 2;

    % Sum all the Pb within the region
    Pb_sum = sum(sum(Pb));
    
    %Method 2 using the center of BBox
    u_x_m = mean(XCoord);
    u_y_m = mean(YCoord);
    
    sigma_x_m = eta*sum(sum(Pb.*(F_XCoord-u_x_m)))/Pb_sum;
    sigma_y_m = eta*sum(sum(Pb.*(F_YCoord-u_y_m)))/Pb_sum;

    G_m = exp(-((F_XCoord-u_x_m).^2./(2*sigma_x_m^2) + (F_YCoord-u_y_m).^2./(2*sigma_y_m^2)));
    ux{i} = u_x_m;
    uy{i} = u_y_m;
    Gm{i} = G_m;
end