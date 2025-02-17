clc;
clear;
m=15;
n=3;
a= cell(m,n);

for i =1:m/3
    for j=1:m/3
      if (i==j)
        a{i,j}=eye(3);
      else
          a{i,j}=zeros(3,3);
      end
    end
end    
        
blk= cell2mat(blkdiag(a))