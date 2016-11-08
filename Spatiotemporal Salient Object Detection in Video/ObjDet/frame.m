str='C:\Documents and Settings\NUS\Desktop\2D3DVideo\2DVideoAndFrames\tree1.avi';
info=aviinfo(str);
fum=info.NumFrames;
for i=1:fum
mov=aviread(str,i);
I=mov.cdata;
    if i < 10
    imwrite(I,strcat('C:\Documents and Settings\NUS\Desktop\2D3DVideo\2DVideoAndFrames\frame1_000',int2str(i),'.bmp'),'bmp');
    elseif i < 100
    imwrite(I,strcat('C:\Documents and Settings\NUS\Desktop\2D3DVideo\2DVideoAndFrames\frame1_00',int2str(i),'.bmp'),'bmp');
    elseif i < 1000
    imwrite(I,strcat('C:\Documents and Settings\NUS\Desktop\2D3DVideo\2DVideoAndFrames\frame1_0',int2str(i),'.bmp'),'bmp');
    elseif i > 999
    imwrite(I,strcat('C:\Documents and Settings\NUS\Desktop\2D3DVideo\2DVideoAndFrames\frame1_',int2str(i),'.bmp'),'bmp');
    end;
end;