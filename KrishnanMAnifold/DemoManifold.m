clc;
clear all;

addpath(genpath(pwd))
load buddha.mat





M=length(cloud);
dim=3;

for i = 1 : M
   transT = rand(dim,1);
   orthoM = compose_rotation(randn(1),randn(1),randn(1));
   cloud_rand{i,1} = cloud{i,1};
   cloud_rand{i,2} =  cloud{i,2}; %orthoM*cloud{i,2} + transT * ones(1,length(cloud{i,2}));%+noise(itrNoise)*randn(dim,length(cloud{i,2}));
end








CorresMatrix=cell(length(cloud_rand),length(cloud_rand));

pts=cell(1,length(cloud_rand));
MInit=zeros(4,4,length(cloud_rand));
for i=1:length(cloud_rand)
MInit(:,:,i)=eye(4);
pts{1,i}=cloud_rand{i,2};
end

for i=1:length(cloud_rand)
    
     for j=1:length(cloud_rand)
     
         if(i~=j)   
                 [~,IA,IB]=intersect(cloud_rand{i,2}',cloud_rand{j,2}','rows');
                 CorresMatrix{i,j}= horzcat(IA,IB);
         end        
     
     end
end


Rmat1=[];
Rmat2=[];
for i=1:length(cloud_rand)
    for j=1:length(cloud_rand)
   
        if ~isempty(CorresMatrix{i,j})
           Rmat1(end+1)=i;
           Rmat2(end+1)=j;
           
        else
            continue;
        end
        
        
    end
end
ijMap=horzcat(Rmat1',Rmat2');


corrMap= cell(length(ijMap),2);%cell(1,length(ijMap));

for i=1:length(ijMap)

corrMap{i}(:,1)=CorresMatrix{ijMap(i,1),ijMap(i,2)}(:,1);
corrMap{i}(:,2)=CorresMatrix{ijMap(i,1),ijMap(i,2)}(:,2);
end

[R,T]=MVSE3EstimateL2Krishnan(MInit,pts,ijMap,corrMap,30,.001);
rot=R;


patches=cell(M,1);
for i=1:M
patches{i,1}=rot(:,((i-1)*3+1):(i*3))*cloud_rand{i,2}+T(:,i)* ones(1,length(cloud_rand{i,2}));
end

b=patches{1,1};
a=patches{2,1};

c=horzcat(patches{:});




