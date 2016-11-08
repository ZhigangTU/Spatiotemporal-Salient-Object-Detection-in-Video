function [frameIOMap, frameObjectMap, frameObjectBoxComb, frameEnergy] = SpTemSaliencyObject(flowdata,magflow,framedata,height,width,nframe,outPath)

currentPath = cd;

options.valScale = 60;
options.alpha = 0.05;
options.color_size = 5;
%% Print status messages on screen
options.vocal = true;
options.regnum =500;
options.m = 20;
options.gradLambda = 1;

videofolder = '001';
options.infolder = fullfile('E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images', 'Diving-Side', videofolder);   %JHMDB\images\brush_hair; data\input\

% The folder where all the outputs will be stored.
options.outfolder = fullfile(currentPath, 'data', 'output\RPCA\Diving-Side', videofolder);
if( ~exist(options.outfolder, 'dir' ) )
    mkdir(options.outfolder )
end
if( ~exist( fullfile( options.outfolder, 'energy'), 'dir' ) )
    mkdir(fullfile( options.outfolder, 'energy'));
end
if( ~exist( fullfile( options.outfolder, 'saliency'), 'dir' ) )
    mkdir(fullfile( options.outfolder, 'saliency'));
end

% % Read all frames in memory
% [data.frames,data.names,height,width,nframe]= readAllFrames(options);

% Load optical flow and frames (or compute if file is not found)
data.frames = framedata;
data.flow = flowdata;

% if(isempty(data.flow ))
%     data.flow = computeOpticalFlow( options, data.frames );
% end

% Load superpixels (or compute if not found)
data.superpixels = loadSuperpixels(options);
if(isempty(data.superpixels))
    data.superpixels = computeSuperpixels(options, data.frames);
end
[superpixels, nodeFrameId, bounds, labels] = makeSuperpixelIndexUnique(data.superpixels);
[colours, centres, t] = getSuperpixelStats(data.frames(1:nframe-1), superpixels, double(labels));

valLAB = [];
for index = 1:nframe-1
    valLAB = [valLAB; data.superpixels{index}.Sup1, data.superpixels{index}.Sup2, data.superpixels{index}.Sup3];     
end
RegionSal = [];
frameEnergy = cell(nframe-1, 1);
frameObject = cell(nframe-1, 1);
frameObjectMap = cell(nframe-1, 1);
frameIOMap = cell(nframe-1, 1);
frameObjectBoxComb = cell(nframe, 1);
frameIOMapBox = cell(nframe, 1);

foregroundArea = 0;

for index = 1:nframe-1        
    frame = data.frames{index};  
    flow  = data.flow{index};
%     frameName = data.names{index};    
    nLabel = max(data.superpixels{index}.Label(:));
    Label = data.superpixels{index}.Label;
    framex = reshape(frame,height*width,1,3);
    Label = reshape(Label,height*width,1);       
    frameVal = colours(bounds(index):bounds(index+1)-1,:);
    framex = uint8(reshape(superpixel2pixel(double(data.superpixels{index}.Label),double(frameVal)),height ,width,3));       
    framex = imfilter(framex,fspecial('average',3),'same','replicate');

%     %% Method 1:Original Consistent Video Saliency Using Local Gradient Flow Optimization
%     % Appearnce Boundary --> Probability of the magnitude of the input image(Not logical labels). Static color RGB boundary-->edge_canny 
%     G = edge_detect(framex); 
%     % Motion magnitude
%     magnitude = magflow{index};
%     if index>1
%         mask = imdilate((frameEnergy{index-1}>0.3),strel('diamond',30))+0.3;
%         mask(mask(:)>1)=1;
%         magnitude = magnitude.*mask;
%         G = G.*mask;
%     end
%     % Probability of the gradient of the flow magnitude
%     gradFlowBoundary = 1 - exp( -options.gradLambda * magnitude );                
%     if (max(magnitude(:))<10)
%         gradFlowBoundary = gradFlowBoundary + 0.01;
%     end
%     % spatio-temporal gradient (RGB Appearance + Flow)
%     if max(magnitude(:)) > 1
%         G = G.*(gradFlowBoundary);         
%         % Classify the gradient probability into 2 categories
%         large = (gradFlowBoundary > 0.55); % Normally 0.55; 0.6
%         G(large) = gradFlowBoundary(large);
%     end
   

    %% Method 2:Appearnce Boundary --> Probability of the magnitude of the input image(Not logical labels). Static appearance boundary --> edge_canny 
    G = 0.1*edge_detect(framex); 
    % Motion Boundary
    magnitude = magflow{index};   
    if index > 1
        mask = imdilate((frameEnergy{index-1}>0.3),strel('diamond',30))+0.3;
        mask(mask(:)>1) = 1;
        magnitude = magnitude.*mask;     % magnitude*0.3, if the frameEnergy of the pixel <= 0.3
        G = G.*mask;                     % G*0.3, if the frameEnergy of the pixel <= 0.3
    end
    % gradBoundary + rotBoundary; + G --> (Whole object) 
    mode = 3;
    if (max(magnitude(:)) >= 1)
        [MoBoundary, gradBoundary] = getProbabilityEdge(flow, magnitude, mode);
        % Probability of the Spatio-temporal gradient (Appearnce Boundary (G) + Motion Boundary (gradBoundary))
        G = 10*G.*(gradBoundary);        % MoBoundary = gradBoundary*rotBoundary (original-->G.*(gradBoundary))              
        % Classify the gradient probability into 2 categories
        large = (gradBoundary > 0.55);   % Normally 0.55; 0.6
        G(large) = gradBoundary(large);
    else
        for Ind = index-1:-1:1
            magnitude = magflow{Ind};
            if max(magnitude(:)) >= 1.0
                flow = data.flow{Ind};
                [MoBoundary, gradBoundary] = getProbabilityEdge(flow, magnitude, mode);
                % Probability of the Spatio-temporal gradient (Appearnce Boundary (G) + Motion Boundary (gradBoundary))
                G = 10*G.*(gradBoundary);       % MoBoundary = gradBoundary*rotBoundary (original-->G.*(gradBoundary))              
                % Classify the gradient probability into 2 categories
                large = (gradBoundary > 0.55);  % Normally 0.6; 0.55
                G(large) = gradBoundary(large);
                break;
            end           
        end      
    end
    frameObject{index} = G; 
    
%     %% Refinement according to inside-outside map
%     Maxmag  = max(magnitude(:))
%     Meanmag = mean(magnitude(:))
%     if max(magnitude(:))>5 && mean(magnitude(:))>0.05
%         T1 = 0.55;
%         T2 = 0.75;
%     elseif (max(magnitude(:))<=5 && max(magnitude(:))>1) && (mean(magnitude(:))<=0.05 && mean(magnitude(:))>0.01)
%         T1 = 0.25;
%         T2 = 0.45;
%     else
%         T1 = 0.10;
%         T2 = 0.25;
%     end
    
%     if (max(magnitude(:)) >= 1)
        GEdge = (G>0.01);   %**G>0.055; 0.006
%     else
%         GEdge = (G>0.55);
%     end 
    % Clean the noisy pixels of the detected Edges 
    GEdgeC = EdgeClean(GEdge, height, width, 0);
    % Performing the InOutMaps -- imdilate(strel('diamond',4)); JHMDB-->4; Sports-->6
    GEdgeMap = getInOutMapsEdge(flow, GEdgeC, 4, outPath, index, 1);  % 4-->insideVotes; 1-->Dilate Lab;  
    frameObjectMap{index} = GEdgeMap;    
    frameObjectBoxComb = BoundingBoxLab(GEdgeMap, index, nframe, frame, frameObjectBoxComb, outPath);  
    
    %% Refine the edge if the object box is large
    TargetsCoord = frameObjectBoxComb{index};
    Row = size(TargetsCoord,1);
    NumLarg = 0;
    for iRow = 1:Row
        iRowRat = TargetsCoord(iRow,6);
        if iRowRat>1/3
            NumLarg = NumLarg+1;
        end
    end
    if NumLarg > 0
        % Refinement according to inside-outside map
%         if (max(magnitude(:)) >= 1) 
            GEdge = (G>0.025);   %**G>0.075
%         else
%             GEdge = (G>0.70);
%         end   
        % Clean the noisy pixels of the detected Edges 
        GEdgeC = EdgeClean(GEdge, height, width, 1);    
        % Performing the InOutMaps
        GEdgeMap = getInOutMapsEdge(flow, GEdgeC, 4, outPath, index, 1);  %4-->insideVotes; 1-->Dilate Lab
        frameObjectMap{index} = GEdgeMap;    
        frameObjectBoxComb = BoundingBoxLab(GEdgeMap, index, nframe, frame, frameObjectBoxComb, outPath);   
    end
        
    GEdgePath = sprintf('%s/GEdgeObject%05d.jpg',outPath, index);    
    imwrite(GEdge, GEdgePath);
    GEdgeMapPath = sprintf('%s/GEdgeObjectC%05d.jpg',outPath, index);    
    imwrite(GEdgeMap, GEdgeMapPath);
%     figure; imshow(GEdge)
%     figure; imshow(GEdgeC)
    figure; imshow(GEdgeMap)
    
    
    %% Method 2 -- Compute the inside-outside map --> Mt  (Motion salient part)
    InOutMap = getInOutMaps(flow, magnitude, outPath, index);  % data.inMaps; Mt = data.inMaps;
    frameIOMap{index} = InOutMap;
    if sum(InOutMap(:))>0
        frameIOMapBox = BoundingBoxLab(InOutMap, index, nframe, frame, frameIOMapBox, outPath, 2);
    else
        frameIOMapBox{index} = [index round(0.4*size(frame,1)) round(0.4*size(frame,2)) round(0.6*size(frame,1)) round(0.6*size(frame,2)) 0];  
    end
    figure; imshow(InOutMap)
    
    
    %% Method 3 -- Learning to Detect Motion Boundaries
    boundaries_Color     = detect_motionboundaries(uint8(frame));      % not accurate
    boundaries_ColorFlow = detect_motionboundaries(uint8(frame), flow);
%     InOutColorMap     = getInOutMapsLDMB(boundaries_Color, outPath, index);
%     InOutColorFlowMap = getInOutMapsLDMB(boundaries_ColorFlow, outPath, index);
%     figure; imshow(boundaries_Color)
%     figure; imshow(boundaries_ColorFlow)
%     figure; imshow(InOutColorMap)
%     figure; imshow(InOutColorFlowMap)
    
    
    %% saliency via gradient flow
    [V_Energy1,H_Energy1,V_Energy2,H_Energy2] = energy_map(double(framex),double(G));

    if index ==1
        Energy = min(min(min(H_Energy1,V_Energy1),H_Energy2),V_Energy2);       
    else
        mask = int32(imdilate((Energy>0.2),strel('diamond',20)));
        mask = ~mask;
        Energymap = (Energy<0.05).*mask; 
        Energymap = ~Energymap;
        Energy = Energy*0.3+(Energymap.*min(min(min(H_Energy1,V_Energy1),H_Energy2),V_Energy2))*0.7;%considering saliency of prior frame
    end

    Energy = Energy/max(Energy(:));         
    L{1} = uint32(data.superpixels{index}.Label);
    S{1} = repmat(Energy,[1 3]);
    [ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, double(nLabel));
    R = double(R(:,1));
    [sR,indexR] = sort(R);
    t = sum(sR(end-9:end))/10;
    R = (R-min(R))/(t-min(R));
    R(R>1)=1;
    RegionSal = [RegionSal;R];
    Energy = reshape(superpixel2pixel(double(data.superpixels{index}.Label),double(R)),height ,width);     
    frameEnergy{index} = Energy;
    
    foregroundArea = foregroundArea + sum(sum(frameEnergy{index}>0.6)); 
    
    EnergyPath = sprintf('%s/EdgeEnergyObject%05d.jpg',outPath, index);    
    imwrite(Energy, EnergyPath);
end

ObjectsCoord = sprintf('%s/frameObjectBoxsRef.mat',outPath);
save(ObjectsCoord,'frameObjectBoxComb');