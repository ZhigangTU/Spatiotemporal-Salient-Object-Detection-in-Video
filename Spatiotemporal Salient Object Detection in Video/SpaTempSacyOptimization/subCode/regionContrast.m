function Consal = regionContrast( valCOL,Label,centres)
    map = zeros(size(Label));
    valSUM = unique( valCOL(:,3)+ 1000*valCOL(:,2)+1000000*valCOL(:,1));
    valSUMall = valCOL(:,3)+ 1000*valCOL(:,2)+1000000*valCOL(:,1);
    colorcentres = centres;
    Consal = ones(size(centres,1),1);
    
    for index =1:size(valSUM,1)
        centre = centres(valSUMall==valSUM(index),:);
        num = size(centre,1);
        centre = sum(centre,1)/num;
        colorcentres(valSUMall==valSUM(index),:) = repmat(centre,num,1);
    end
    Consal = sqrt(sum((centres-colorcentres).^2,2));
    Consal(Consal < 0.5) = 0.5;
    t = max(Consal(:))*0.5;
    Consal = (Consal-0.5)/(t-0.5);
    Consal(Consal >1) = 1;
end