function [xreturn] = nl_conjgrad_fluor_test(A, b, x,niter,varargin)

% TESTING 
%1 only l1 norm on PFOB image/
%2: l2 norm on sum of 
visualizationoption=1; 

realI=varargin{1};
lambda0=varargin{2};
n1=varargin{3};
n2=varargin{4};
visualizationoption=varargin{5};


rr = @(I) reshape(I,[n1,n2]);
T= MakeWaveletOp(x,n1,n2);
% T=opDirac(n1*n2)
x=T*x;
A=A*(T'); 
mask=abs(b)>0; 

if visualizationoption
figure(100);clf 
subplot(221); imshow(rr(T'*x),[]); drawnow; end

lambda=0; %test:
grad=-gradient(A,b,x,T,lambda,mask);
s=grad;
t0=1;
bb=0.7;
for i=1:niter
    if i>3; lambda=lambda0; end
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
            disp(['iter:',num2str(i)]); xreturn=T'*x; 
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
%   beta=0 ;% Newton??
    
    % update conjugate direction:
    sold=s;
    s=grad+beta*sold;
    
    if visualizationoption
    figure(100); 
    subplot(211);
    hold on 
    imshow(rr(abs(T'*x)),[]); colormap jet
    text(5,5,num2str(i),'Color','white')
    hold off
    drawnow;

    subplot(212); 
    hold on; plot(i,f1,'k.')
    plot(i,fl1,'ro')
    plot(i,fl2,'g*')
    hold off; drawnow;
    pause(1);
    title('obj (bl); l1 (red); l2(green)')
    
    fprintf('iter: %d |  obj: %d | l2: %d | l1: %d \n',(i),f1,fl2,fl1)

    end
        

end

xreturn=T'*x; %return back to image space!

end

function grad=gradient(A,b,x,T,lambda,mask) % for image domain at least; 

gradl1=gradientl1(A,b,x,T);
grad=2*A'*(gather(mask).*(A*x-gather(b)))+lambda*gradl1; 

end

function [objective,fl1,objectivel2]=Calcobjective(A,b,s,T,lambda,mask)
obj=(mask.*(A*s-b));
objectivel2=obj(:)'*obj(:);

if 0
    fprintf('showing l2 norm: %f \n',objectivel2)
    As=A*s; 
    figure(201); plot(abs(As(6.5e4:6.6e4))); hold on; plot(abs(b(6.5e4:6.6e4))); hold off;
    pause(.1);
end

fl1=lambda*sum(objectivel1(A,b,s,T));

objective=sum(objectivel2+fl1);
end


function objective=objectivel1(A,b,s,T)
if 0; s(1:length(s)/2)=0; end;
objective=((s).*conj(s)+eps).^(1/2); 
end

function gradl1=gradientl1(A,b,x,T)
if 0; x(1:length(x)/2)=0; end;
gradl1 = x.*(x.*conj(x)+eps).^(-1/2);
end

function W=MakeWaveletOp(x,n1,n2)
W=opWavelet2(n1,n2,'Haar',4,4,0);
end
