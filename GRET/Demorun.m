%
% main file computing the GRET-SDP, GRET-SPEC and GRET-ADMM 
% input data is the correspondences and the 3D point clouds 

clc;
clear all;

addpath(genpath(pwd));
load('bunmain.mat')
pcl_lib =1 ; % if pcl library is present then 1
             % else 0;

M=length(cloud);    % 'cloud' is the cell containing patches of 
                    %  point clouds
                    %  length of cloud gives total number of point clouds

dim=3;              % dimension ofpoint cloud 
N=length(bun);      % total number of distinct points in the cloud

j=0;
numItr=2;          % total number of iteration to run the code


% add random rotation and translation to the cloud
% add controlled random noise to the point cloud

noise=0:0.00025:0.009;


mainCloud=horzcat(cloud{:,2});
 
savepcd('bun.pcd',bun');
savepcd('bigbun.pcd',mainCloud);



for itr=1:numItr
     j=0;
for itrNoise=1:length(noise)       %                   noise=0;%.0010:0.009
    fprintf('iteration %d  noise level %d noise %f\n',itr,j,noise(itrNoise));
    j=j+1;
    resultMat(j,1,itr)=noise(itrNoise);
    tic
    cloud_rand={length(cloud),2}; % add random rotation and translation to the cloud

    
   trans_rand=cell(M);
   rot_rand =cell(M);
 for i = 1 : M
 
   transT = noise(itrNoise)*randn(dim,1); 
   trans_rand{i}=transT;
   orthoM = compose_rotation(rand(1),rand(1),rand(1));
   rot_rand{i}=orthoM;
   cloud_rand{i,1} = cloud{i,1};
   cloud_rand{i,2} = orthoM*cloud{i,2} + transT * ones(1,length(cloud{i,2}))+noise(itrNoise)*randn(dim,length(cloud{i,2}));
end


% Crete the L,D,B matrix for noisy point cloud
[L,D,B]=createLDBMatrix2(cloud_rand,N);

% calculate the final position of each point along with rotations and translation of each
% patch

[X2,rADMM,tADMM]=GRET_ADMM_vin(cloud_rand,N,L,D,B);
savepcd('temp1.pcd',X2);

%%
[X3,rSpec,tSpec]=GRET_SPEC_vin(cloud_rand,N,L,D,B);
savepcd('temp2.pcd',X3);


%% Krishnan



pts=cell(1,length(cloud_rand));
MInit=zeros(4,4,length(cloud_rand));
for i=1:length(cloud_rand)
MInit(:,:,i)=eye(4);
pts{1,i}=cloud_rand{i,2};
end

[R,TKris]=MVSE3EstimateL2Krishnan(MInit,pts,ijMap,corrMap,30,.001);
rot=R;

rKrish=cell(M,1);
patches=cell(M,1);
for i=1:M
rKrish{i}=rot(:,((i-1)*3+1):(i*3));    

patches{i}=rKrish{i}*cloud_rand{i,2}+TKris(:,i)* ones(1,length(cloud_rand{i,2}));
end

X4=horzcat(patches{:});










 %% Error Calc
 
 
















 savepcd('X3.pcd',X3);
 unix('./svdr X3.pcd bun.pcd');
 trmat1=importtrmat('example.txt')  ;
 rot=trmat1(1:3,1:3);
 trans= trmat1(1:3,4);
 tmp3=X3+trans*ones(1,length(X3));
 resultMat(j,2,itr)= errorANE((rot*(tmp3))',bun');
 
  resultMat(j,2,itr)=0;

 



 
 
 savepcd('X4.pcd',X4);
 unix('./svdr X4.pcd bigbun.pcd');
 trmat1=importtrmat('example.txt');
 rot=trmat1(1:3,1:3);
 trans= trmat1(1:3,4);
 tmp3=X4+trans*ones(1,length(X4'));
 
resultMat(j,4,itr)=  errorANE((rot*(tmp3))',mainCloud);

 
 
 resultMat(j,2,itr)=0;
 resultMat(j,4,itr)=0;
 
 
 savepcd('X2.pcd',X2);
 unix('./svdr X2.pcd bun.pcd');
 trmat1=importtrmat('example.txt');
 rot=trmat1(1:3,1:3);
 trans= trmat1(1:3,4);
 size(X2);
 tmp3=X2+trans*ones(1,length(X2));
resultMat(j,3,itr)=errorANE((rot*(tmp3))',bun');
 
 
 
 
 
 
 
 
 


fprintf('%.6f %.6f %.6f \n',resultMat(j,2,itr),resultMat(j,3,itr),resultMat(j,4,itr));


toc
end
end

% save the results along with time and date
savefile=sprintf('Results_%s.mat',datestr(now)); 
save(savefile, 'resultMat');

%% plot the graph

specErr=zeros(length(resultMat),1);
admmErr=zeros(length(resultMat),1);
sdpErr=zeros(length(resultMat),1);


for k=1:length(resultMat)
specErr(k)=sum(resultMat(k,2,:))/numItr;
sdpErr(k)=sum(resultMat(k,4,:))/numItr;
admmErr(k)=sum(resultMat(k,3,:))/numItr;
end

%xax=0:0.0010:0.009; % noise levels added



% plot the semilogy graph 

hold on
semilogy(noise',specErr,'b');
semilogy(noise',admmErr,'r');
semilogy(noise',sdpErr,'m');
grid on;

hold off

%% plot the corresponding 3D-bunny 



if(pcl_lib==1)
    pclviewer(X1);
    waitforbuttonpress;
    pclviewer(X2);
    waitforbuttonpress;
    pclviewer(X3);
else


Xspec=X1';
Xadmm=X2';
Xsdp=X3';

scatter3(Xspec(:,1),Xspec(:,2),Xspec(:,3));
waitforbuttonpress
scatter3(Xadmm(:,1),Xadmm(:,2),Xadmm(:,3));
waitforbuttonpress
scatter3(Xsdp(:,1),Xsdp(:,2),Xsdp(:,3));
end
