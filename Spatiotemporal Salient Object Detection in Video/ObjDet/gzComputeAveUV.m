function [aveU,aveV] = gzComputeAveUV(Pframecur, uframecur, vframecur,curH, curW, nSize)
    [nH, nW] = size(Pframecur);
    nhStart = (curH - 1)*nSize + 1;
    nwStart = (curW - 1)*nSize + 1;
    OneBlockpixelMask = Pframecur(nhStart : min(nhStart + nSize -1, nH), nwStart : min(nwStart + nSize -1, nW));
    OneBlockpixelMask = OneBlockpixelMask>0.6;

    OneBlocku = uframecur(nhStart : min(nhStart + nSize -1, nH), nwStart : min(nwStart + nSize -1, nW));
    OneBlockv = vframecur(nhStart : min(nhStart + nSize -1, nH), nwStart : min(nwStart + nSize -1, nW));
    OneBlocku = OneBlocku.*OneBlockpixelMask;
    OneBlockv = OneBlockv.*OneBlockpixelMask;

    num  = sum(OneBlockpixelMask(:));
    aveU = sum(OneBlocku(:))/num;
    aveV = sum(OneBlockv(:))/num;

end


% nhStart : min(nhStart + nSize -1, nH), nwStart : min(nwStart + nSize -1, nW)