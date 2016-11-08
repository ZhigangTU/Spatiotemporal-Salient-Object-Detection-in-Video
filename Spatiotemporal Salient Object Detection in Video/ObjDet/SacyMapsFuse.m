function [SpaTempConSacyMap,RefinConSacyMapsPCA,RefinConSacyMapsFOS] = SacyMapsFuse(optFgConSacyMapsPCA,optFgConSacyMapsFOS,TMapsPCA,TMapsFOS)

if TMapsPCA < 0.25
    TAMapsPCA = TMapsPCA*1.25;
else
    TAMapsPCA = TMapsPCA;
end
% fprintf('\nTMapsPCA=%3.3f TAMapsPCA=%3.3f FraNum=%d\n', TMapsPCA, TAMapsPCA, j);
% outObjectPath = sprintf('%s/FgConSacyMapsPCASegTAu%05d.jpg',outPath, j);
% optFgConSacyMapsPCATLab = (optFgConSacyMapsPCA>=TAMapsPCA);
% imwrite(optFgConSacyMapsPCATLab, outObjectPath);
    
if TMapsFOS>0.50
    TAMapsFOS = TMapsFOS*0.80;
elseif TMapsFOS<0.25
    TAMapsFOS = TMapsFOS*1.25;  %(1/0.8)
else
    TAMapsFOS = TMapsFOS;
end
% fprintf('\nTMapsFOS=%3.3f TAMapsFOS=%3.3f FraNum=%d\n', TMapsFOS, TAMapsFOS, j);
% outObjectPath = sprintf('%s/FgConSacyMapsFOSSegTAu%05d.jpg', outPath, j); 
% optFgConSacyMapsFOSTLab = (optFgConSacyMapsFOS>=TAMapsFOS);
% imwrite(optFgConSacyMapsFOSTLab, outObjectPath);
    
Ra = 0.50; % [1.00,0.75,0.50]    
optFgConSacyMapsPCATLab = (optFgConSacyMapsPCA>TAMapsPCA*Ra); % [1.00,0.75,0.50]
optFgConSacyMapsFOSTLab = (optFgConSacyMapsFOS>TAMapsFOS*Ra);

OvlapPixs = optFgConSacyMapsPCATLab.*optFgConSacyMapsFOSTLab;
OvlapPixs = OvlapPixs>0;
if sum(OvlapPixs(:)) > 0
    % Initialization
    RefinConSacyMapsPCA = optFgConSacyMapsPCA;
    RefinConSacyMapsFOS = optFgConSacyMapsFOS;
    % Refinement
    [mP, nP] = size(optFgConSacyMapsPCATLab);
    TargetPAC = optFgConSacyMapsPCATLab>0.6;
    TargetFOS = optFgConSacyMapsFOSTLab>0.6;
    [LPCA,NumPCA]=bwlabel(TargetPAC,8);
    [LFOS,NumFOS]=bwlabel(TargetFOS,8);
    for iTargFOS = 1:NumFOS
        [rtemp,ctemp] = find(LFOS==iTargFOS);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);
        minH = max(1,minH);
        minW = max(1,minW);
        maxH = min(mP,maxH);
        maxW = min(nP,maxW);
        % TargetOBox = [minH minW maxH maxW];
        box1 = [minW,minH,(maxW-minW),(maxH-minH)];   %[X,Y,WIDTH,HEIGHT]
        area1 = (maxW-minW)*(maxH-minH);
        for iTargPCA = 1:NumPCA
            [Prtemp,Pctemp] = find(LPCA==iTargPCA);
            PminH=min(Prtemp);
            PminW=min(Pctemp);
            PmaxH=max(Prtemp);
            PmaxW=max(Pctemp);
            PminH = max(1,PminH);
            PminW = max(1,PminW);
            PmaxH = min(mP,PmaxH);
            PmaxW = min(nP,PmaxW);
            box2 = [PminW,PminH,(PmaxW-PminW),(PmaxH-PminH)];
            area2 = (PmaxW-PminW)*(PmaxH-PminH);
            iou = inters_union(box1,box2);
            if iou > 0
                if iou > 0.75
                    if area1<area2
                        RefinConSacyMapsFOS(PminH:PmaxH,PminW:PmaxW) = 2*optFgConSacyMapsFOS(PminH:PmaxH,PminW:PmaxW);
                    else
                        RefinConSacyMapsFOS(minH:maxH,minW:maxW) = 2*optFgConSacyMapsFOS(minH:maxH,minW:maxW);
                    end
                else
                    RefinConSacyMapsFOS(minH:maxH,minW:maxW) = 2*optFgConSacyMapsFOS(minH:maxH,minW:maxW);
                end
            end
        end
    end
    for iTargPCA = 1:NumPCA
        [rtemp,ctemp] = find(LPCA==iTargPCA);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);
        minH = max(1,minH);
        minW = max(1,minW);
        maxH = min(mP,maxH);
        maxW = min(nP,maxW);
        box1 = [minW,minH,(maxW-minW),(maxH-minH)];   %[X,Y,WIDTH,HEIGHT]
        area1 = (maxW-minW)*(maxH-minH);
        for iTargFOS = 1:NumFOS
            [Frtemp,Fctemp] = find(LFOS==iTargFOS);
            FminH=min(Frtemp);
            FminW=min(Fctemp);
            FmaxH=max(Frtemp);
            FmaxW=max(Fctemp);
            FminH = max(1,FminH);
            FminW = max(1,FminW);
            FmaxH = min(mP,FmaxH);
            FmaxW = min(nP,FmaxW);
            box2 = [FminW,FminH,(FmaxW-FminW),(FmaxH-FminH)];
            area2 = (FmaxW-FminW)*(FmaxH-FminH);
            iou = inters_union(box1,box2);
            if iou > 0
                if iou > 0.75
                    if area1<area2
                        RefinConSacyMapsPCA(FminH:FmaxH,FminW:FmaxW) = 2*optFgConSacyMapsPCA(FminH:FmaxH,FminW:FmaxW);
                    else
                        RefinConSacyMapsPCA(minH:maxH,minW:maxW) = 2*optFgConSacyMapsPCA(minH:maxH,minW:maxW);
                    end
                else
                    if area1 > area2
                        RefinConSacyMapsPCA(FminH:FmaxH,FminW:FmaxW) = 2*optFgConSacyMapsPCA(FminH:FmaxH,FminW:FmaxW);
                    else
                        RefinConSacyMapsPCA(minH:maxH,minW:maxW) = 2*optFgConSacyMapsPCA(minH:maxH,minW:maxW);
                    end
                end
            end
        end
    end
%     RefinConSacyMapsFOS(OvlapPixs) = 2*RefinConSacyMapsFOS(OvlapPixs);
%     RefinConSacyMapsPCA(OvlapPixs) = 2*RefinConSacyMapsPCA(OvlapPixs);
else
    RefinConSacyMapsFOS = optFgConSacyMapsFOS;
    RefinConSacyMapsPCA = optFgConSacyMapsPCA;
end

% Normalization
maxSacyMapsFOS = max(RefinConSacyMapsFOS(:));
maxSacyMapsPCA = max(RefinConSacyMapsPCA(:));
RefinConSacyMapsFOS = RefinConSacyMapsFOS/maxSacyMapsFOS;
RefinConSacyMapsPCA = RefinConSacyMapsPCA/maxSacyMapsPCA;


% Refinement the 2th time 
TReMapsFOS = graythresh(RefinConSacyMapsFOS);
TReMapsPCA = graythresh(RefinConSacyMapsPCA);
RefinConSacyMapsFOSTLab = (RefinConSacyMapsFOS>=TReMapsFOS);
RefinConSacyMapsPCATLab = (RefinConSacyMapsPCA>=TReMapsPCA);
OvlapPixsRe = RefinConSacyMapsFOSTLab.*RefinConSacyMapsPCATLab;
OvlapPixsRe = OvlapPixsRe>0;
SpaTempComLab = zeros(size(RefinConSacyMapsFOSTLab));
if sum(OvlapPixsRe(:))>0
    [mP, nP] = size(RefinConSacyMapsFOSTLab);
    TargetPAC = optFgConSacyMapsPCATLab>0.6;
    TargetFOS = optFgConSacyMapsFOSTLab>0.6;
    [LPCA,NumPCA]=bwlabel(TargetPAC,8);
    [LFOS,NumFOS]=bwlabel(TargetFOS,8);
    for iTargFOS = 1:NumFOS
        [rtemp,ctemp] = find(LFOS==iTargFOS);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);
        minH = max(1,minH);
        minW = max(1,minW);
        maxH = min(mP,maxH);
        maxW = min(nP,maxW);
        box1 = [minW,minH,(maxW-minW),(maxH-minH)];   %[X,Y,WIDTH,HEIGHT]
        area1 = (maxW-minW)*(maxH-minH);
        for iTargPCA = 1:NumPCA
            [Prtemp,Pctemp] = find(LPCA==iTargPCA);
            PminH=min(Prtemp);
            PminW=min(Pctemp);
            PmaxH=max(Prtemp);
            PmaxW=max(Pctemp);
            PminH = max(1,PminH);
            PminW = max(1,PminW);
            PmaxH = min(mP,PmaxH);
            PmaxW = min(nP,PmaxW);
            box2 = [PminW,PminH,(PmaxW-PminW),(PmaxH-PminH)];
            area2 = (PmaxW-PminW)*(PmaxH-PminH);
            iou = inters_union(box1,box2);
            if iou > 0.75
%                 if area1 < area2
%                     SpaTempComLab(FminH:FmaxH,FminW:FmaxW) = 1;
%                 else
%                     SpaTempComLab(minH:maxH,minW:maxW) = 1;
%                 end
                CMinH = min(minH,PminH);CMinW = min(minW,PminW);CMaxH = max(maxH,PmaxH);CMaxW = max(maxW,PmaxW);
                SpaTempComLab(CMinH:CMaxH,CMinW:CMaxW) = 1;
            end
        end
    end
    % Identify the overlapped regions
    SpaTempComLab(OvlapPixsRe) = 1;
    % Combing Appearance-Motion SacyMap
    SpaTempOvrLab = (SpaTempComLab>0);
    SpaTempConSacyMap(SpaTempOvrLab) = max(RefinConSacyMapsPCA(SpaTempOvrLab), RefinConSacyMapsFOS(SpaTempOvrLab));
    SpaTempConSacyMap(~SpaTempOvrLab) = RefinConSacyMapsPCA(~SpaTempOvrLab).*RefinConSacyMapsFOS(~SpaTempOvrLab);
    SpaTempConSacyMap = reshape(SpaTempConSacyMap, size(SpaTempOvrLab));
else
    % Combing Appearance-Motion SacyMap
    SpaTempConSacyMap = RefinConSacyMapsPCA.*RefinConSacyMapsFOS;
end

% Normalize the Spatio-temperal Saliency Maps
maxSpaTempSacyMap = max(SpaTempConSacyMap(:));
SpaTempConSacyMap = SpaTempConSacyMap/maxSpaTempSacyMap;

TSpaTempMaps = graythresh(SpaTempConSacyMap);
TSpaTempMaps = max(TSpaTempMaps*(Ra/4),0.1);

SpaTempConSacyMap(SpaTempConSacyMap<TSpaTempMaps) = eps;


% CMinH = min(minH,FminH);CMinW = min(minW,FminW);CMaxH = max(maxH,FmaxH);CMaxW = max(maxW,FmaxW);
% RefinConSacyMapsPCA(CMinH:CMaxH,CMinW:CMaxW) = 2*optFgConSacyMapsPCA(CMinH:CMaxH,CMinW:CMaxW);
