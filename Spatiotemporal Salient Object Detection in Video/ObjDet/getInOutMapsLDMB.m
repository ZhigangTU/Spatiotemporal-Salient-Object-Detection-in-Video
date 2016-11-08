function output = getInOutMapsLDMB(boundaryMap, outPath, frameIndx)   
    [height, width, ~] = size(boundaryMap);
    
    % Motion boundaries touching the edges will be cut!
    sideCut = false(height, width);
    sideCut(1: 20, :) = true;
    sideCut(end - 20: end, :) = true;
    sideCut(:, 1: 20) = true;
    sideCut(:, end - 20: end) = true;

    insideVotes = getInPoints(boundaryMap, sideCut, false);     % original(false); (boundaryMap / flowGraBoundary)

    % Logical Labels
    outObject = insideVotes > 4;   % >4 = 5

    % Delect and clean noisy pixels
    [L,numtarget] = bwlabel(outObject,8);  
    stats = regionprops(L,'Area');  % Compute the size of each connected region
    area = cat(1,stats.Area);     
    index = find(area > 5^2);      % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    outObjectC = ismember(L,index(:));

    output = outObjectC;

%     outObjectPath = sprintf('%s/VotesObject%05d.jpg',outPath, frameIndx);    
%     imwrite(outObject, outObjectPath);
%     outObjectCPath = sprintf('%s/VotesObjectC%05d.jpg',outPath, frameIndx);    
%     imwrite(outObjectC, outObjectCPath); 
    
end

