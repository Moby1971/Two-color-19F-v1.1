% draft: implement irwls in 3d recon
% we have a recon RCG

function k = BlindDeconvolutionOptimization(RCG,N,linear_recon,nx,Spectrum1,rr,varargin)
%varargin(1)= directions
%varargin(2)= separation option
%varargin(2)= visualization option

if nargin>6
    directions=varargin{1};
    nacq=length(directions);
    if nargin>7
        sep=varargin{2}; %separate spectrum per acquistion
    else
        sep=1;
    end
    
    if nargin>8 
        visualization=varargin{3};
    else
        visualization=1;
    end
    
else
    directions={'FH','RL'};
    nacq=2;
    visualization=1;
    sep=0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFOB=RCG(1:N);
PFCE=RCG(N+1:2*N);

x_orig=PFOB;        % assumed original image (most recent PFOB recon);
Y=linear_recon-repmat(PFCE,[nacq,1]);     %convolved measurement (lin_recon minus PFCE);

%%%%%DEFINE LINEAR OPERATOR X
x0=rr(x_orig); % original image

%% make circ shift matrices
for jj=1:nacq
    X1=[];
    if strcmp(directions{jj},'LR')
        for nn=1:nx
            xtemp=circshift(x0,[0 nn-1]);
            X1=cat(2,X1,xtemp(:));
        end
        
    elseif strcmp(directions{jj},'RL')
        for nn=1:nx
            xtemp=circshift(x0,[0 -nn]);
            X1=cat(2,X1,xtemp(:));
        end
        
    elseif strcmp(directions{jj},'FH')
        for nn=1:nx
            xtemp=circshift(x0,[nn-1 0]);
            X1=cat(2,X1,xtemp(:));
        end
        
    elseif strcmp(directions{jj},'HF')
        for nn=1:nx
            xtemp=circshift(x0,[-nn 0]);
            X1=cat(2,X1,xtemp(:));
        end
    else
        fprintf('Direction not found! \n')
    end
    
    MC{jj}=X1;
end



%%

Z=zeros(size(X1));
X=[];

% X Z Z Z ; Z X Z Z; Z Z X Z; Z Z Z X %block diagonal matrix
if sep  %separate spectrum per acquisition
    for ii=1:nacq
        R=[repmat(Z,[1 (ii-1)]),MC{ii},repmat(Z,[1 nacq-(ii)])]; %row
        X=[X;R];
    end
else %optimize one spectrum for all acquisitions
    for ii=1:nacq
        R=[MC{ii}]; %row
        X=[X;R];
    end
end

%%

%%% solve
opts.pcg_tol = 1e-8;
opts.pcg_its = 2;
opts.lambda=1e6;

%     validregion=validregion+circshift(validregion,-1)+circshift(validregion,1)

if sep
    k_init=repmat(Spectrum1,[1,nacq]); % USING PREVIOUS ESTIMATE FOR PSF
    % k_init=zeros(1,nacq*nx); % NOT USING PREVIOUS ESTIMATE FOR PSF
else
    k_init=Spectrum1;
    %     k_init=zeros(1,nx);
end

%validregion=k_init>0; %T

for iter=1:5 %outer iterations
    fprintf('Blind Deconv outer iter:%4.2i \n',iter)
    k=IRWLS_BlindDeconvolution(k_init.', X, Y(:),opts);
    
    %normalize and threshold;
    k((k) < 0) = 0;
    k((k) < 0.01*max(abs(k(:)))) = 0; %remove noise 

%     k(~validregion)=0;
    %         k(k<0.02*max(k(:)))=0;
    sumk = sum((k(:))); 
    k = k ./ sumk;
    if sep; k=nacq*k; end;
    
    k_init=k.';
    
    if visualization
    figure(70);clf
    hold on;
    if ~sep
        plot(abs(k),'+-'); plot(abs(Spectrum1),'o');
    else
        plot(abs(k),'+-'); plot(abs(repmat(Spectrum1,[1 nacq])),'o');
    end
    hold off; drawnow; title('original and optimized kernel ')
    pause(0.01)
    end
end


if sep==0
   k=repmat(k,[1 nacq]); %repeat output;  
end

end