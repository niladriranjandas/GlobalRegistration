function trans = recover_trans(T,E,rots)


    n = length(rots); m =length(E);
    I = [1:m,1:m]'; J=zeros(2*m,1);
    
%     S=[-1*ones(m,1);ones(m,1)];
   
    for i=1:m
        J(i) = I(i);
        J(i+m) = J(i);
        T(i,:)=rots{E(i)}*T(i,:)';
    end
  
    A=zeros(m,n);
    
     for i=1:length(E)
       A(i,E(i,1))=-1;
       A(i,E(i,2))=1;
         
    end
% %      disp(I);
% %      disp(J);
% %     A=sparse(I,J,S,m,n);
%    % A=zeros(m,n);
%    
%      A(E+1,1)=1;
%     T(:,E+1)=0;
% %     
% disp('T');
% disp(T);
   A=sparse(A);
    t1=A\T(:,1);
    t2=A\T(:,2);
    t3=A\T(:,3);
    trans=[t1,t2,t3]; %trans should be 9X3 matrix
%     
%     trans=T;
end

      
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   