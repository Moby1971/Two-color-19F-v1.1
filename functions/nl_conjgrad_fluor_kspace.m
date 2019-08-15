function [xreturn] = nl_conjgrad_fluor_kspace(A, b, x,niter,varargin)


realI=varargin{1};
realflag=1;
lambda=varargin{2}
n1=varargin{3}
n2=varargin{4};
kspace_input=varargin{5};
FS=varargin{6}; %block diagonal fourier op

rr = @(I) reshape(I,[n1,n2])
FBlock=opBlockDiag(FS,FS)
T= MakeWaveletOp(x,n1,n1)
T=opBlockDiag(T,T)*FBlock';

figure(100);clf 
% x=zeros(size(A'*b)); 
subplot(221); imshow(rr(FBlock'*x)); drawnow;
mask=abs(b)>0; 

grad=-gradient(A,b,x,T,lambda,mask);
s=grad;

t0=1;
bb=0.7
for i=1:niter
    t=t0;
    % perform line search
    
    lsiter=0;
    alpha=1;
    f0=Calcobjective(A,b,x,T,lambda,mask);
    f1=Calcobjective(A,b,x+alpha*s,T,lambda,mask);

    
    while ((f1) > (f0 - alpha*t*abs(s(:)'*grad(:))))  % change this 
        [f1,fl1,fl2]=Calcobjective(A,b,x+alpha*s,T,lambda,mask);
        alpha=alpha.*t*bb;
        if lsiter>100; disp('line search failed, convergence reached?');
            
            if realflag==1
                mse=sum((x-realI).^2);
            else
                mse=0;
            end
            
            disp(['iter:',num2str(i),' MSE: ',num2str(mse)])
            
            return; end
        
        lsiter=lsiter+1;
    end
    
    if lsiter > 2
		t0 = t0 * bb;
	end 
	
	if lsiter<1
		t0 = t0 / bb;
	end
    disp(['lsiter=',num2str(lsiter),' alpha=',num2str(alpha)])
    % update the position
    xold=x; 
    x=xold+alpha*s;
    
    % calculate steepest direction
    gradold=grad;
    grad=-gradient(A,b,x,T,lambda,mask);
    
    % calculate beta: TO DO 
    beta=(grad'*(grad-gradold))/(gradold'*gradold+eps); %POLAK-RIBIERE 
    %beta=0 ;% Newton??
    
    % update conjugate direction:
    sold=s;
    s=grad+beta*sold;
    
    figure(100); 
    subplot(221);
    hold on 
    imshow(rr(abs(FBlock'*x)),[]); drawnow;
    text(5,5,num2str(i),'Color','white')
    hold off
    drawnow;
    
    if realflag==1
    mse=sum((x-realI).^2);
    else
        mse=0;
    end
    subplot(222); 
    hold on; plot(i,f1,'k.')
    plot(i,fl1,'ro')
    plot(i,fl2,'g*')
    hold off; drawnow;
    pause(1);
    title('obj (bl); l1 (red); l2(green)')
    
    
    fprintf('iter: %d | MSE: %d | obj: %d | l2: %d | l1: %d \n',(i),(mse),f1,fl2,fl1)

end

xreturn=x; %return back to image space!

end

function grad=gradient(A,b,x,T,lambda,mask) % for image domain at least; 

gradl1=gradientl1(A,b,x,T);
grad=2*A'*(mask.*(A*x-b))+lambda*gradl1; 
figure(1); clf ; hold on; plot(abs((A*x)-b),'k-'); plot(abs(b),'b.'); plot(abs(A*x),'r+'); drawnow; pause(1)
end

function [objective,fl1,objectivel2]=Calcobjective(A,b,s,T,lambda,mask)
obj=(mask.*(A*s-b));
objectivel2=obj(:)'*obj(:);
fl1=lambda*sum(objectivel1(A,b,s,T));

objective=sum(objectivel2+fl1);
end


function objective=objectivel1(A,b,s,T)
s=T*s; 
objective=((s).*conj(s)+eps).^(1/2); 

figure(2); clf; imshow(reshape(abs(s(:)),128,256),[0 1]);
end

function gradl1=gradientl1(A,b,x,T)
x=T*x; 
gradl1 = x.*(x.*conj(x)+eps).^(-1/2);
gradl1=T'*gradl1;
end

function W=MakeWaveletOp(x,n1,n2)

W=opWavelet2(n1,n2,'Daubechies',4,4,0);

end
