
% clc;
% clear all;
  calc_rot_corr;
 
 



[Mmat,TransMat,Emat,Dmat,wVect]=constructRot(pairs_rot,wtMat);

E=Emat;
w=wVect;
numIter=3;



[rots,trans]=TSR_global_align(Mmat,TransMat,E,w,numIter);
tmp=[0 0 0 1];
transForm=cell(length(rots));
for i=1:length(rots)
transForm{i}=horzcat(rots{i},trans(i,:)');
transForm{i}=vertcat(transForm{i},tmp);
end



homeDir ='home/vinith';
dataTemp='dataN';
dirMain='~/GbReg';
dirName = '~/GbReg/data';%# folder path
dataDir   ='GbReg/data';
binDirName='GbReg/trimesh2/bin.Linux64';
binDir='trimesh2/bin.Linux64';

for i=1:length(file)

file1=sprintf('%s/%s/%s',dirName,dataTemp,file{i});
file2=sprintf('~/%s',binDirName);    

trp=sprintf('cp %s %s',file1,file2);
unix(trp);
[~,fname,ext]=fileparts(file{i,1});
filename=sprintf('%s.xf',fname);
dlmwrite(filename,transForm{i},'delimiter',' ','precision',3)
file3=sprintf('%s/%s',dirMain,filename);
trp=sprintf('cp  %s %s',file3,file2);
unix(trp);

createNply = sprintf('(cd ../%s && ./mesh_xform %s %s %s_main.ply 0 0 )',binDirName,file{i,1},filename,fname);

unix(createNply);






end



cleanxf = sprintf('(cd ../%s && rm -r *.xf )',binDirName);
unix(cleanxf);
cleanply = sprintf('(cd ../%s && rm -r *_N.ply )',binDirName);
unix(cleanply);
cleanxf = sprintf('rm -r *.xf ');
unix(cleanxf);


%cd(binDir);

mply= dir( fullfile(binDir,'*main.ply') );
emt=' ';
for i=1:length(mply)
   sng=mply(i).name;
    emt=sprintf('%s %s',emt,sng);
    
end


dispply = sprintf('(cd ../%s && ./mesh_view %s )',binDirName,emt);
unix(dispply);
dispply = sprintf('(cd ../%s && ./mesh_cat %s -o combineEIG.ply )',binDirName,emt);
unix(dispply);

cd(binDir);

dispply = sprintf('cp %s ~/GbReg','combineEIG.ply');
unix(dispply);
unix('rm -r *.ply');
unix('rm -r *.pair');
unix('rm -r *.xf');

cd ..;
cd ..;




disp(pwd);

%% 
Mmat_=SDPmat(Mmat);
[rots,trans]=TSR_global_align(Mmat_,TransMat,E,w,numIter);
tmp=[0 0 0 1];
transForm=cell(length(rots));
for i=1:length(rots)
transForm{i}=horzcat(rots{i},trans(i,:)');
transForm{i}=vertcat(transForm{i},tmp);
end



homeDir ='home/vinith';
dataTemp='dataN';
dirMain='~/GbReg';
dirName = '~/GbReg/data';%# folder path
dataDir   ='GbReg/data';
binDirName='GbReg/trimesh2/bin.Linux64';
binDir='trimesh2/bin.Linux64';

for i=1:length(file)

file1=sprintf('%s/%s/%s',dirName,dataTemp,file{i});
file2=sprintf('~/%s',binDirName);    

trp=sprintf('cp %s %s',file1,file2);
unix(trp);
[~,fname,ext]=fileparts(file{i,1});
filename=sprintf('%s.xf',fname);
dlmwrite(filename,transForm{i},'delimiter',' ','precision',3)
file3=sprintf('%s/%s',dirMain,filename);
trp=sprintf('cp  %s %s',file3,file2);
unix(trp);

createNply = sprintf('(cd ../%s && ./mesh_xform %s %s %s_main.ply 0 0 )',binDirName,file{i,1},filename,fname);

unix(createNply);






end



cleanxf = sprintf('(cd ../%s && rm -r *.xf )',binDirName);
unix(cleanxf);
cleanply = sprintf('(cd ../%s && rm -r *_N.ply )',binDirName);
unix(cleanply);
cleanxf = sprintf('rm -r *.xf ');
unix(cleanxf);


%cd(binDir);

mply= dir( fullfile(binDir,'*main.ply') );
emt=' ';
for i=1:length(mply)
   sng=mply(i).name;
    emt=sprintf('%s %s',emt,sng);
    
end


dispply = sprintf('(cd ../%s && ./mesh_view %s )',binDirName,emt);
unix(dispply);
dispply = sprintf('(cd ../%s && ./mesh_cat %s -o combineSDP.ply )',binDirName,emt);
unix(dispply);
disp(pwd);

cd(binDir);

dispply = sprintf('cp %s ~/GbReg','combineSDP.ply');
unix(dispply);
unix('rm -r *.ply');
unix('rm -r *.pair');
unix('rm -r *.xf');

cd ..
cd ..
% 
% 
% DispICPerr = sprintf('./mesh_align recon.ply combineICP.ply');
% DispEIGerr = sprintf('./mesh_align recon.ply combineEIG.ply');
% DispSDPerr = sprintf('./mesh_align recon.ply combineSDP.ply');
% 
% 
% [~,ICPerr]=unix(DispICPerr); fprintf('The error using sequentical ICP : %f',str2double(ICPerr));
% [~,EIGerr]=unix(DispEIGerr); fprintf('The error using sequentical EIG : %f',str2double(EIGerr));
% [~,SDPerr]=unix(DispSDPerr); fprintf('The error using sequentical SDP : %f',str2double(SDPerr));
% 

disp('The End');
