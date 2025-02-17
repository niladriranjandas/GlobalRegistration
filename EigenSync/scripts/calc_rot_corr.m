

clc;clear all;


homeDir ='home/vinith';
dataTemp='dataN';
dirName = '~/GbReg/data';%# folder path
dataDir   ='GbReg/data';
OrgRecon  ='GbReg/recon.ply';
binDirName='GbReg/trimesh2/bin.Linux64';
binDir='trimesh2/bin.Linux64';
files = dir( fullfile(dirName,'*.ply') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
pairs= cell(length(files),length(files));
pairs_rot= cell(length(files),length(files));
wtMat=zeros(length(files),length(files));
Patch=cell(length(files),1);
noise=0.0100;
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

clc;
temp=sprintf('%s/%s',dirName,dataTemp);
file = dir( fullfile(temp,'*_N.ply') );
file={file.name}';

% disp(pwd);
for i=1:length(file)

file1=sprintf('%s/%s/%s',dirName,dataTemp,file{i});
file2=sprintf('~/%s',binDirName);    
% disp(file1);
% disp(file2);
trp=sprintf('cp %s %s',file1,file2);
unix(trp);
end

cd(binDir);
% check if no ply already exist in the folder

for i=1:length(file)
         extractVtx = sprintf('./mesh_print %s',file{i});
         [stat,cmd1]=unix(extractVtx); % if error here ...change permission to make ./mesh_print executable
         Patch{i}=importvertex('hello.vtx');
         delete('hello.vtx');
         
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




cd ..;
cd ..;

PatchMat = Patch{1,1};
for i = 2:length(Patch)

PatchMat = vertcat(PatchMat,Patch{i,1});
end

n1=length(wtMat);
wtMat = wtMat(1:n1-1,1:n1-1);
pairs_rot = pairs_rot(1:n1-1,1:n1-1);
pairs = pairs(1:n1-1,1:n1-1);

disp(pwd);














%% Seq ICP begins
for i=1:length(file)

file1=sprintf('%s/%s/%s',dirName,dataTemp,file{i});
file2=sprintf('~/%s',binDirName);    
% disp(file1);
% disp(file2);
trp=sprintf('cp %s %s',file1,file2);
unix(trp);
end
cd(binDir);

MemGraph=cell(length(file),length(file));
disp(pwd);
% check if no ply already exist in the folder
fileLeft=1:length(file);
clc;
for i=1:length(fileLeft)
   i1=1;
   for j=1:length(fileLeft)
       
       wt=ones(1,length(fileLeft));
       disp(wt);
         if (i~=j)
             
             createICP = sprintf('./mesh_align_combine %s %s',file{i1},file{j});
             [~,name1,~]=fileparts(file{i1});[~,name2,~]=fileparts(file{j});
            [~,cmdout]=unix(createICP); dist = str2double(cmdout); 
            wt(j)=dist;
             disp(dist);
             if dist~=0
             pairfile = sprintf('%s__%s.pair',name1,name2);
             
             delete(pairfile);
             pairfile = sprintf('%s_%s.xf',name1,name2);
             delete(pairfile);
          
                
             end
       
             
         end  
         
   end
   
   
             [q,w]=min(wt);
             createICP = sprintf('./mesh_align_combine %s %s',file{i1},file{w});
             [~,name1,~]=fileparts(file{i1});[~,name2,~]=fileparts(file{j});
             [status,cmdout]=unix(createICP); dist = str2double(cmdout); 
             pairfile = sprintf('%s__%s.pair',name1,name2);
             MemGraph{i,j}=importpairfile(pairfile);
             delete(pairfile);
             createNply = sprintf('./mesh_xform %s.ply %s_%s.xf %s.ply %f %f',name2,name1,name2,name2,noise,FlagNoise);
             unix(createNply);disp(createNply);
             createCAT = sprintf('./mesh_cat %s.ply %s.ply -o %s.ply',name1,name2,name1);
             unix(createCAT);
             fileLeft(w)=[];
%              
             
end

 dispply = sprintf('./mesh_view %s ',file{1});
 unix(dispply);

movefile(file{1},'combineICP.ply');


unix('rm -r *.xf');
% unix('rm -r *.pair');















for i=2:length(file)
delete(file{i});
end




unix('rm -r *.xf');
% unix('rm -r *.pair');


dispply = sprintf('cp %s ~/GbReg','combineICP.ply');
unix(dispply);




cd ..;
cd ..;

disp(pwd);

clearvars -except pairs pairs_rot wtMat Patch Patch PatchMat edgeMat file MemGraph;
