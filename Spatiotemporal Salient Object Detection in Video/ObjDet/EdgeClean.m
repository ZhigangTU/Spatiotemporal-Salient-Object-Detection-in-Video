function ReEdge = EdgeClean(GEdge, height, width, VoteLab)

hsz = 2;
% if(~exist('VoteLab', 'var') || isempty(VoteLab))
%     PixelT = hsz*2;
% else    
%     PixelT = hsz*2+1;
% end
PixelT = hsz*2+VoteLab;

sz = [height, width];
[indx_row, indx_col] = find(GEdge ==1);      
pad_im = padarray(GEdge, hsz*[1 1], 'symmetric', 'both');            
[H, W] = size(pad_im);

ReEdge = GEdge;
% Divide into several groups for memory reasons ~70,000 causes out of memory
Indx_Row = indx_row;
Indx_Col = indx_col;
N        = length(Indx_Row); % number of elements to process
n        = 1e4;              % number of elements per batch
nB       = ceil(N/n);

for ib = 1:nB;
    istart = (ib-1)*n + 1;
    iend   = min(ib*n, N);
    indx_row = Indx_Row(istart:iend);
    indx_col = Indx_Col(istart:iend);    

    [C, R] = meshgrid(-hsz:hsz, -hsz:hsz);
    nindx = R + C*H;    
    cindx = indx_row +hsz  + (indx_col+hsz-1)*H;
    
    pad_indx = repmat(nindx(:), [1 length(indx_row)]) + repmat(cindx(:)', [(hsz*2+1)^2, 1] );
    sumEdge = sum(pad_im(pad_indx),1);
    sumEdgeL = (sumEdge>PixelT);
    
    indx = sub2ind(sz, indx_row, indx_col);
    ReEdge(indx) = sumEdgeL;
end



% [indx_row, indx_col] = find(GEdge ==1);
% Npix = length(indx_row);
% half_win = 1;
% for ii=1:Npix
%     IndxR = indx_row(ii);
%     IndxC = indx_col(ii);
%     minrow = max(IndxR-half_win,1);
%     maxrow = min(IndxR+half_win,height);
%     mincol = max(IndxC-half_win,1);
%     maxcol = min(IndxC+half_win,width);
%     sumEg = 0;
%     for r=minrow:maxrow
%         for c=mincol:maxcol
%             sumEg = sumEg+GEdge(r,c);
%         end
%     end
%     if sumEg<(2*half_win+1)
%        GEdge(indx_row,indx_col)=0; 
%     end
% end
% ReEdge = (GEdge>0.5);