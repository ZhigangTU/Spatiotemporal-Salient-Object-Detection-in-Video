function output = getInOutMapsEdge(flow, boundaryMap, TVotes, outPath, index, DiaLab)

    if nargin < 6
        DiaLab = false;
    end
               
    [height, width, ~] = size(flow);
    
    % Motion boundaries touching the edges will be cut!
    sideCut = false(height, width);
    sideCut(1: 20, :) = true;
    sideCut(end - 20: end, :) = true;
    sideCut(:, 1: 20) = true;
    sideCut(:, end - 20: end) = true;

    insideVotes = getInPoints(boundaryMap, sideCut, false);     % original(false); (boundaryMap / flowGraBoundary)
 
%     % Refinement the insideVotes based on flow
%     edgeT = 0.5;    % 0.5; 0.25; 0.1
%     if(getFrameQuality(insideVotes > 4) < 0.2)  % a pixel with 5 or more rays intersecting the boundaries an odd number of times is deemed inside
%         boundaryMap = calibrateProbabilityEdge(flow, 0.71, edgeT);
%         insideVotes = getInPoints(boundaryMap, sideCut, false); % original(false)
%     end

    % Logical Object Labels
    outObject = insideVotes > TVotes;   % >4 = 5
    % Clean the noisy pixels of the detected Edges
    outObject = EdgeClean(outObject, height, width, 0);
    
%     GEdgePath = sprintf('%s/GEdgeObjectVotes%05d.bmp',outPath, index);    
%     imwrite(outObject, GEdgePath); 
    
    if DiaLab
        % Delect and clean noisy pixels
        outObjectD = imdilate(outObject,strel('diamond',6));  %[JHMDB-->(4,12^2); Sports-->(6,20^2)]
        [L,numtarget]=bwlabel(outObjectD,8);      
        stats = regionprops(L,'Area');  % Compute the size of each connected region
        area = cat(1,stats.Area);     
        index = find(area > 20^2);      % Finding the index of the region that its size larger than a threshold (e.g.10*10); [JHMDB-->12^2,Sports-->20^2]
        % Obtaining the indexed regions which eliminate small noisy regions 
        outObjectC = ismember(L,index(:));

        output = outObjectC;
        
    else
        output = outObject;
    end

end



%     outputObject = sprintf('%s/Object%05d.bmp',outPath, frameN);    
%     imwrite(outObject, outputObject);
%     outputObjectD = sprintf('%s/ObjectD%05d.bmp',outPath, frameN);    
%     imwrite(outObjectD, outputObjectD);
%     outputObjectC = sprintf('%s/ObjectC%05d.bmp',outPath, frameN);    
%     imwrite(outObjectC, outputObjectC); 