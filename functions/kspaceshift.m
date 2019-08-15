function K= kspaceshift(K,d) 
r1=linspace(-d(1),d(1),size(K,1));
r2=linspace(-d(2),d(2),size(K,2)); r2=r2.';

K=bsxfun(@times,K,exp(-1i*r1)) ;
K=bsxfun(@times,K,exp(-1i*r2)) ;
end

