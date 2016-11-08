function [wCtr, thresh] = SaliencyObjectnessTu(TargetLabels,H,W,idxImg,pixelList,adjcMatrix,colDistM,clipVal)

N = length(pixelList);
% Set the number of proposals to consider
% For fast mode, set nProposals to 200. For more accurate results but at a slower rate, set nProposals to 1000.
% nProposals = 1000;   

% % Reads Objectness Proposals
% A = dlmread([ 'BingBoxes/' bbox_filename '.txt']);
% X = A(2:end,:);

% X = bbox_filename;
% MAT = zeros(H,W);
% spMAP = zeros(N,1);
% 
% for i = 1:nProposals
%     hwin = X(i,4)-X(i,2)+1;
%     wwin = X(i,5)-X(i,3)+1;
%     wMat = gausswin(hwin,2)*gausswin(wwin,2)'; 
%     for x = X(i,3):X(i,5)
%         for y = X(i,2):X(i,4)
%             spMAP(idxImg(y,x)) = spMAP(idxImg(y,x)) + X(i,1)*wMat(y-X(i,3)+1,x-X(i,2)+1);
%        end
%     end   
% end

for PTresh = 0.65:-0.05:0.60   %(0.75:-0.05:0.50)
    isFG = [];
    for i = 1:N
        MAT = zeros(H,W);
        MAT(pixelList{i,1}) = 1;     % MAT=(idxImg==i), MAT is a super-pixel which contains the pixels in it
        fspNumP = sum(MAT(:));
        spFG = MAT.*TargetLabels;    % Finding the index of the super-pixel locates on the foreground object candidates 
        isFGspi = sum(spFG(:))>PTresh*fspNumP; % For one super-pixel, if it has more than PTresh% percent pixels locates on candidates, it is considered as foreground super-pixel.
        isFG = [isFG;isFGspi];
    end
    if sum(isFG)>=100
        break;
    end 
end
% FinalPTresh = PTresh

geoDistMatrix = CalGeoDist(adjcMatrix, colDistM, clipVal);

FGscore = zeros(N,1);
BGscore = zeros(N,1);

for i = 1:N,
    for j = 1:N,
        if isFG(j) == 1,
            FGscore(i) = FGscore(i) + geoDistMatrix(i,j);
        else
            BGscore(i) = BGscore(i) + geoDistMatrix(i,j);
        end
    end
end

wCtr = zeros(N,1);
MAT  = zeros(H,W);
for i = 1:N
    MAT(pixelList{i,1}) = BGscore(i)/(FGscore(i)+ eps);
    wCtr(i) = BGscore(i)/(FGscore(i)+ eps);    % weighted contrast-->wCtr
end

wCtr = (wCtr - min(wCtr)) / (max(wCtr) - min(wCtr) + eps);
IwCtr = wCtr;
thresh = graythresh(wCtr);      % automatic threshold

% Method 1 (Original)
% if thresh>0
%     wCtr(wCtr < thresh) = 0;  % The original method to threshold the Objectness Map to obtain the salient object
% else
%     wCtr(isFG==0) = 0;
% end
    
% % Method 2 (Our)
% if thresh>0  
% %     IwCtr = wCtr;
% %     wCtr(wCtr < min(thresh+0.25,0.80)) = 0;
%     wCtr(wCtr < min(thresh*2.0,0.80)) = 0;
% %     wCtr(isFG>0) = IwCtr(isFG>0);
% else
%     wCtr(isFG==0) = 0;
% end

% % Method 3 (Our)
% wCtr(isFG==0) = 0;

% Method 4 (Our)
if thresh > 0
    mT = min(max(median(wCtr(isFG>0))+0.1,0.80),0.90); %+0.10
    mTs = max(median(wCtr)+0.1,0.10); % mTs = max(median(wCtr(isFG>0))-0.1,0.10);
    wCtrs = wCtr;   
    wCtr(wCtr < mT) = 0; 
    wCtr(isFG>0) = IwCtr(isFG>0);     % IwCtr(isFG>0)||wCtr(wCtr < mT)
    wCtrs(wCtrs < mTs) = 0;
    for j = 1:N
        if wCtr(j)>0 && wCtrs(j)==0   % Minus wCtrs(wCtrs < mTs)
            wCtr(j) = 0;
        end
    end
else
    wCtr(isFG==0) = 0;
end

AutGrayThresh = thresh;






% if thresh == 0
%     LwCtr  = length(wCtr);
%     LBeCtr = length(BeCtr);
%     if LwCtr>LBeCtr
%         wCtr = zeros(LwCtr);
%         wCtr(1:LBeCtr) = BeCtr;
%     else
%         wCtr = BeCtr(1:LwCtr);
%     end
% end

%     Obj1 = zeros(H,W);
%     Obj2 = zeros(H,W);
%     for i = 1:N
%         if wCtr(i)>0
%             Obj1(pixelList{i,1}) = 1;
%         end
%     end
%     IwCtr(isFG==0) = 0;
%     for i = 1:N
%         if IwCtr(i)>0
%             Obj2(pixelList{i,1}) = 1;
%         end
%     end