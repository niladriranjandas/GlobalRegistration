
function [Mmat,TransMat,Emat,Dmat,wVect] = constructRot(RotTrans,wt)

temp = length(RotTrans);

Mmat =cell(temp,temp);
trans=cell(temp,temp);



for i=1:(length(RotTrans))
    for j=1:(length(RotTrans))
      
        if i==j
        Mmat{i,j}=eye(3,3);
        trans{i,j} =  zeros(3,1);
        else
            disp(i);disp(j);
            if (iscell(RotTrans{i,j}(1:3,1:3))==1)
                    Mmat{i,j}=cell2mat(RotTrans{i,j}(1:3,1:3));
            else
                     Mmat{i,j}=RotTrans{i,j}(1:3,1:3);
            end
            
            if( iscell(RotTrans{i,j}(1:3,4))==1)
              trans{i,j} = cell2mat(RotTrans{i,j}(1:3,4));
            else
              trans{i,j} = RotTrans{i,j}(1:3,4);
            end
        end
        
        
    Mmat{i,j}=wt(i,j)*Mmat{i,j};    
    end
end

Mmat=cell2mat(Mmat);



numel=nnz(wt);

TransMat=zeros(numel,3);

[row,col,val]=find(wt);

Emat=[row col];
wVect=val;
for i= 1:length(Emat)
   
TransMat(i,:) =(trans{Emat(i,1),Emat(i,2)});
end

dmat =sum(wt,2);

Dmat=cell(length(dmat),length(dmat));

for i = 1:(length(dmat))
    for j = 1:(length(dmat))
    
     if i==j 
         Dmat{i,j}=dmat(i)*eye(3);
     else 
         Dmat{i,j}=zeros(3,3);
     end
    end
end

Dmat= cell2mat(Dmat);















end