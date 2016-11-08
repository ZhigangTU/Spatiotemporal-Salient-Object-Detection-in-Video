MATT = zeros(H,W);
for j = 1:N
    if wCtr(j)>0
        MATT(pixelList{j,1}) = 1;
    end
end
figure; imshow(MATT)
