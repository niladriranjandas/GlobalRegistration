
clc;clear all;


homeDir ='home/vinith';
dataTemp='dataN';
dirName = '~/GbReg/data';%# folder path
dataDir   ='GbReg/data';
binDirName='GbReg/trimesh2/bin.Linux64';
binDir='trimesh2/bin.Linux64';
files = dir( fullfile(dirName,'*.ply') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
pairs= cell(length(files),length(files));
pairs_rot= cell(length(files),length(files));
wtMat=zeros(length(files),length(files));
Patch=cell(length(files),1);
noise=0.004;
FlagNoise=0;


% -------------------------------------------------------------------
% get the conf file
%
%extract data convert the xf and then transform the original patch 
% and store the transform patch in dataN
%--------------------------------------------------------------------

conFile= dir( fullfile(dirName,'*.conf') );
conFile ={conFile.name}';
confName =fullfile(dirName,conFile{1});
RotTrans = importConf(confName);
for i=1:length(files)-1
[~,fname,ext]=fileparts(RotTrans{i,2});disp(fname);
% createXf=sprintf('(cd ../%s && ./xf -trans %f %f %f -q %f %f %f %f > %s.xf)',binDirName,RotTrans{i,3},RotTrans{i,4},RotTrans{i,5},RotTrans{i,6},RotTrans{i,7},RotTrans{i,8},RotTrans{i,9},fname);
createXf=sprintf('(cd ../%s && ./xf -trans %f %f %f  -v %f %f %f %f > %s.xf)',binDirName,RotTrans{i,3},RotTrans{i,4},RotTrans{i,5},RotTrans{i,6},RotTrans{i,7},RotTrans{i,8},RotTrans{i,9},fname);
unix(createXf);
fnameTemp = sprintf('%s/%s.xf',binDir,fname);
copyfile(fnameTemp,dirName);
createNply = sprintf('(cd ../%s && ./mesh_xform %s/%s.ply %s/%s.xf %s/%s/%s_N.ply %f %f)',binDirName,dirName,fname,dirName,fname,dirName,dataTemp,fname,noise,FlagNoise);
% disp(createNply);
[a,b{i}]=unix(createNply);
delete(fnameTemp);
end
% disp(fnameTemp);
% -------------------------------------------------------------------
%% 

clc;
temp=sprintf('%s/%s',dirName,dataTemp);
file = dir( fullfile(temp,'*_N.ply') );
file={file.name}';

for i=1:length(file)

file1=sprintf('%s/%s/%s',dirName,dataTemp,file{i});
file2=sprintf('~/%s',binDirName);    
% disp(file1);
% disp(file2);
trp=sprintf('cp %s %s',file1,file2);
unix(trp);
end

cd(binDir);


%% 
Fs=1:length(file);
wt=1:length(file);
for i=1:length(file)
%          extractVtx = sprintf('./mesh_print %s',file{i});
%          [stat,cmd1]=unix(extractVtx); % if error here ...change permission to make ./mesh_print executable
%          Patch{i}=importvertex('hello.vtx');
%          delete('hello.vtx');
         
   for j=1:length(file)
         if (i~=j)
             
             createICP = sprintf('./mesh_align_combine %s %s',file{i},file{j});
             [~,name1,~]=fileparts(file{i});[~,name2,~]=fileparts(file{j});
             [status,cmdout]=unix(createICP); dist = str2double(cmdout); 
             wtMat(i,j)= dist; %1/(1+dist);
             disp(dist);
             if dist~=0
             pairfile = sprintf('%s__%s.pair',name1,name2);
             pairs{i,j}=importpairfile(pairfile);
             delete(pairfile);
             pairfile = sprintf('%s_%s.xf',name1,name2);
             pairs_rot{i,j}=importRot(pairfile,1,3);
             delete(pairfile);
             else
                 pairs{i,j}=0; pairs_rot{i,j}=zeros(3,4); wtMat(i,j)=0;
             end
         else
              pairs{i,j}=0; pairs_rot{i,j}=zeros(3,4); wtMat(i,j)=0;
         end
    end
end

for i=1:length(file)
delete(file{i});
end


