
function [MSE_R,MSE_T]=MSE(true_R,estim_R,true_T,estim_T)
   
      n=length(true_R);
      Q=estim_R'*true_R/n;
      s=svd(Q);
      MSE_R=6-2*sum(s);
      
      
      tmean =mean(true_T,2);
      emean  = mean(estim_T,2);
      
      for i=1:n
          MSE_T = MSE_T + norm(true_T(:,i)-tmean-estim_T(:,i)+emean);
      end
      MSE_T=MSE_T/n;
end