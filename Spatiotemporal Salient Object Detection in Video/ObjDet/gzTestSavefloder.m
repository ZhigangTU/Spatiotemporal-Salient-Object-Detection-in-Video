function gzTestSavefloder(data)
%%
%path
currentPath = cd;
userName = 'savetest' ;
destDir = fullfile(currentPath,userName) ;
if ~exist(destDir,'dir')
    mkdir(currentPath,userName) ;
end
%%
[m n] = size( data );
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);
%%
for j = 1:5,
  imagedata   = reshape( data(:,j), nH, nW );
  outputNum  = sprintf('test%05d.bmp',j);
  outputName = fullfile(destDir, outputNum);
  imwrite(imagedata, outputName);
end

end
