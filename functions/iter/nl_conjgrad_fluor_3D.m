function [xreturn] = nl_conjgrad_fluor_3D(app, A, b, x, niter, nouter, outeriter, T, lambda, size_im,  varargin)

% TESTING 
%1 only l1 norm on PFOB image/
%2: l2 norm on sum of 

visualizationoption=0;
n1=size_im(1); 
n2=size_im(2); 
n3=size_im(3); 
slice=round(n3/2);
realI=[]; 

if nargin==11
    visualizationoption=varargin{1};
elseif nargin==12
    slice=varargin{2}; 
    visualizationoption=varargin{1};
elseif nargin==13
    realI=varargin{3}; 
    slice=varargin{2}; 
    visualizationoption=varargin{1};
end

rr = @(I) reshape(I,[n1,n2,n3]);
rr1 = @(I) reshape(I,[n1,n2,n3,2]);

if isempty(T)
    T=opDirac(n1*n2*n3)
end

mask=abs(b)>0; 



visualize(x,slice,rr1,visualizationoption)

grad=-gradient(A,b,x,T,lambda,mask);
s=grad;
t0=1;
bb=0.7;

for i=1:niter

    TextMessage(app,strcat("CG iteration ",num2str(i)," ... "));

    t=t0;
    
    % perform line search
    lsiter=0;
    alpha=1;
    f0=Calcobjective(A,b,x,T,lambda,mask);
    f1=Calcobjective(A,b,x+alpha*s,T,lambda,mask);
    
    
    while (gather(f1) > gather(f0 - alpha*t*abs(s(:)'*grad(:))))  % change this
        [f1,fl1,fl2]=Calcobjective(A,b,x+alpha*s,T,lambda,mask);
        alpha=alpha.*t*bb;
        
        if lsiter>100
            
            TextMessage(app,"Line search failed, convergence reached ?");
            TextMessage(app,strcat("Number of iteration = ",num2str(i)));
            
            return;
        end
        
        lsiter=lsiter+1;
    end
    
    if lsiter > 2
        t0 = t0 * bb;
    end
    
    if lsiter<1
        t0 = t0 / bb;
    end
    
    TextMessage(app,strcat("lsiter = ",num2str(lsiter),"   alpha = ",num2str(alpha)));
    
    % update the position
    xold=x;
    x=xold+alpha*s;
    
    % calculate steepest direction
    gradold=grad;
    grad=-gradient(A,b,x,T,lambda,mask);
    
    % calculate beta: TO DO
    beta=(grad(:)'*(grad(:)-gradold(:)))/(gradold(:)'*gradold(:)+eps); %POLAK-RIBIERE
    %   beta=0 ;% Newton??
    
    % update conjugate direction:
    sold=s;
    s=grad+beta*sold;
    
    
    % Show the fit progress in the GUI
    progress = round(100*((outeriter-1)*niter + i)/(nouter * niter));
    ShowFitProgress(app,progress);
    
    % Shows the intermediate results of the reconstruction in the GUI
    ShowReco(app,rr1(x));
    
    % Display the fitting parameters in the GUI
    TextMessage(app,strcat("Obj = ",num2str(f1,'%10.2e'),"  l2 = ",num2str(fl2,'%10.2e'),"  l1 = ",num2str(fl1,'%10.2e')));
    
    
end

xreturn=x; %return back to image space!

end

%======================================================================

function grad = gradient(A,b,x,T,lambda,mask) % for image domain at least; 

gradl1=gradientl1(A,b,x,T);
grad=2*A'*(mask.*(A*(x)-b))+lambda*gradl1;   

end

%======================================================================

function [objective,fl1,objectivel2]=Calcobjective(A,b,s,T,lambda,mask)

obj=(mask.*(A*s-b));
objectivel2=real(obj(:)'*obj(:));

fl1=lambda*sum((objectivel1(A,b,s,T)));
objective=sum(objectivel2(:)+fl1(:));

end

%======================================================================


function objective=objectivel1(A,b,s,T)

% if 1; s(1:length(s)/2)=0; end;  %% L1 only on PFCE????
s=T*s; %do TV operator 

objective=((s(:)).*conj(s(:))+eps).^(1/2); 

end

%======================================================================

function gradl1=gradientl1(A,b,x,T)

s=size(x); 
x=T*x(:); 

% if 1; x(1:length(x)/2)=0; end; %?????? WHAT IS THIS
gradl1 = x.*(x.*conj(x)+eps).^(-1/2);
gradl1=T'*gradl1;

gradl1=reshape(gradl1,s);

end

%================ Not used for now... 2D wavelet case =================

function W=MakeWaveletOp(x,n1,n2)

W=opWavelet2(n1,n2,'Haar',4,4,0);

end

%======================================================================

function visualize(x,slice,rr1,visualizationoption)

%%% FUNCTION OBSOLETE IN THE APP %%%%

if visualizationoption
    x=rr1(x);
    figure(1001);
    hold on
    imshow(cat(2,squeeze(abs(x(:,:,slice,1))),squeeze(abs(x(:,:,slice,2)))),[]); colormap jet
    hold off
    drawnow;
end

end



