clear all; close all; clc;
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\imagine'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor'))


%%

folder='L:\basic\divi\Projects\cosart\fluor\fluor\recon_realdata\dec17';
protonimages='protonimage_1129.mat'
kspaces='kspaces_Mouse2_48H.mat'


doubleFOV=0;        % 1: double FOV in phase encoding dims/ 0: half FOV in readout dir
zerofilloption=0

fprintf('F-19 reconstruction... \n')
cd(folder)
fprintf('loading kspace data... \n')
load(protonimages)
load(kspaces)

% preprocessing the data
run fix_double_FOV.m

% initialize reconstruction cell
Recon{1,1}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{1,2}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{2,1}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{2,2}=zeros(size(K1,1),size(K1,2),size(K1,3));


% reconstruction options:
dim= 3;             % shared phase_encoding dimension(??) 
phaseremoval=0;     % option to remove phase 
phasecheckerboardremoval=~phaseremoval; 

use_CG=1;           %use CG recon

senselikePFCE=0;
normthreshold=.99;

use_empirical_alpha=1; 
empirical_alpha=[64 101 64 175 64]; % peak heights PFOB
empirical_alpha=[184 502 235 427 274]

visualization=1; 

% parameters 
BW=4.4643e+04;
ppm=282.5685;

run recon_set_params_and_variables.m

%%
%relevant parameters 

lambda=200;
niter=25;
n_outer=3;
IRWLS=0;
opts.pcg_tol = 1e-8;
opts.pcg_its = 100;
opts.lambda=0;


for sl=1:64
tic;   
fprintf('reconstructing slice: %i \n', sl)

k1=extract_slice(ifft_meas(K1,dim),sl,3);
k3=extract_slice(ifft_meas(K3,dim),sl,3);

data=[vec(k1); vec(k3)]; % data matrix 

if phaseremoval
fprintf('removing phase of image data... \n')
data= F2*abs(opInverse(F2)*data); % 
end

if phasecheckerboardremoval
fprintf('removing phase checkerboard image data... \n')
data= (opInverse(F2)*data); % sl
checkerboard1=(-1).^linspace(1,ny,ny).*1i;
checkerboard2=(-1).^linspace(1,nx*2,nx*2).*1i;
data=rr(data);
data=bsxfun(@times,data,checkerboard1.');
data=bsxfun(@times,data,checkerboard2);
data=F2*data(:);
end

k1=data(1:N);
k3=data(N+1:2*N);

% linear recon + visualization 
fprintf('linear (zero-filled) recon... \n')
linear_recon=opInverse(F2)*data;    %linear recon of data 


%%%%%%%%%% visualize linear reconstruction 
if visualization; 
figure(1);subplot(211); imshow(abs(rr(linear_recon)),[]); title('magnitude linear recon') 
figure(1);subplot(212); imshow(angle(rr(linear_recon)),[-pi, pi]); title('phase linear recon');
end

% reshaped magnitude of linear images
l1=(rr1(linear_recon(1:N)));
l3=(rr1(linear_recon((N)+1:2*(N))));


% make convolution operators (include fftshift in operators for simplicity)
Spectrum_vert=Spectrum1.';
Spectrum=(Spectrum1);
A1=opConvolve(nx,ny,Spectrum_vert,[1 1],'cyclic');
A3=opConvolve(nx,ny,Spectrum1,[1 1],'cyclic'); %perhaps 65?
B=opDirac(nx*ny); 


if senselikePFCE
    fprintf('Recon: Sense-like PFCE (MORE TESTING NEEDED!) /n')
    
    l1norm=abs(l1(:));
    l3norm=abs(l3(:));
    sos=sqrt(l1norm.^2+l3norm.^2+1e-6);
    allnorms=cat(3,l1norm./sos,l3norm./sos);
    normstochange=(max(allnorms,[],3)>normthreshold);
    
    if visualization
    figure(22);imshow(rr1(normstochange)); end

    allnorms(normstochange,:,:)=ones(sum(normstochange),1,2)./2;
    sos=sqrt(sum((allnorms.^2),3)+1e-6);%recalculate SoS
    
    % to do : find out which of these options is the best
    % B1=opDiag(allnorms(:,1,1)./sos);
    % B3=opDiag(allnorms(:,1,2)./sos);
    B1=opDiag(2*l1./(l1+l3));
    B3=opDiag(2*l3./(l1+l3));
else
    B1=B;  B3=B;
end

deconv_image1=opInverse(A1)*(opInverse(FS)*data(1:N));
deconv_image3=opInverse(FS*A3)*data((N)+1:2*(N));
Dirac_image1=opInverse(B)*(opInverse(FS)*data(1:N));

%%%%%%%%%% visualize Deconvolution 
if visualization
SpectrumImage1=bsxfun(@times,Spectrum_vert,ones(size(l1))).*max(l1(:));
SpectrumImage2=bsxfun(@times,Spectrum1,ones(size(l1))).*max(l3(:));
figure(5); clf
imshow(abs(cat(2,cat(1,l1,rr1(deconv_image1),SpectrumImage1),cat(1,l3,rr1(deconv_image3),SpectrumImage2))),[]); title('deblurred images')
colormap jet
end

% image-space recon ('first-guess')
M1=[FS*A1,FS*B1];
M3=[FS*A3,FS*B3];
M=[M1;M3];
first_guess=pinv(M)*data;

%%%%%%%%%% visualize CG result 
if visualization
figure(20); imshow(rr2(abs(first_guess)),[]); axis off; title('first guess'); colormap('jet'); end
Recon{2,1}(:,:,sl)=rr1(abs(first_guess(1:N)));
Recon{2,2}(:,:,sl)=rr1(abs(first_guess(N+1:2*N)));

%% iterative optimization 
A1=opConvolve(nx,ny,Spectrum_vert,[1 1],'cyclic');
A3=opConvolve(nx,ny,Spectrum1,[1 1],'cyclic'); %perhaps 65?
M1=[FS*A1,FS*B1];
M3=[FS*A3,FS*B3];
M=[M1;M3];


CG_input=zeros(size(first_guess));
    
if use_CG
    k_init=[Spectrum1,Spectrum1];

    for outer_iter=1:n_outer

    pinvsol=pinv(M)*data;
        
    fprintf('CG recon... \n')
    solutionsize=[2*N,1];

    
    RCG=nl_conjgrad_fluor_test(M,data,CG_input,niter,zeros(solutionsize),lambda,nx,ny*2,[],[]);
    CG_input=RCG;
    Recon{1,1}(:,:,sl)=rr1(abs(RCG(1:N)));
    Recon{1,2}(:,:,sl)=rr1(abs(RCG(N+1:2*N)));
    
    
    if IRWLS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %      OPTIMIZE SPECTRUM ADD HERE 
    PFOB=RCG(1:N);
    PFCE=RCG(N+1:2*N);
%     PFOB=first_guess(1:N);
%     PFCE=first_guess(N+1:2*N);

    x_orig=PFOB;        % assumed original image (most recent PFOB recon);
    Y=linear_recon-cat(1,PFCE,PFCE);     %convolved measurement (lin_recon minus PFCE); 
    
    figure(70); clf; 
    subplot(221); imshow(cat(2,abs(rr1(PFOB)),abs(rr1(PFCE))),[]); title('original PFOB')
    subplot(222); imshow(abs(rr(Y)),[]); title('convolved PFOB measurements')
    subplot(425); plot(abs(k_init)); title('spectrum magnitude')
    subplot(427); plot(angle(k_init)); title('spectrum phase')

    %%%%%DEFINE LINEAR OPERATOR X 
    x0=rr1(x_orig); % original image
    X1=[]; X2=[];
    for nn=1:sqrt(N)
        xtemp=circshift(x0,[0 nn-1]);
        X1=cat(2,X1,xtemp(:));
    end
    
    for nn=1:sqrt(N)
        xtemp=circshift(x0,[nn-1 0]);
        X2=cat(2,X2,xtemp(:));
    end
    
    Z=zeros(size(X1)); 
    X=[X2,Z;Z,X1];
    
    k_init=zeros(1,2*nx); % NOT USING PREVIOUS ESTIMATE FOR PSF
    for iter=1 %outer iterations
        disp(iter)
        k=IRWLS_BlindDeconvolution(k_init.', X, Y(:),opts);
        
        %normalize and threshold;
        k((k) < 0) = 0;
%         k(~validregion)=0;        
%         k(k<0.02*max(k(:)))=0;
        sumk = sum(k(:));
        k = 2*k ./ sumk;
        
        k_init=k.';
        
        figure(70); 
        subpl4=subplot(224); cla(subpl4);
        hold on;
        plot(abs(k),'+-'); plot(abs([Spectrum1,Spectrum1]),'o');
        hold off; drawnow; title('original and optimized kernel ')
    end
    
    
    
    % update convolution and measurement operator 
    
    A1=opConvolve(nx,ny,k(1:nx),[1 1],'cyclic');
    A3=opConvolve(nx,ny,k(nx+1:end).',[1 1],'cyclic'); %perhaps 65?
    M1=[FS*A1,FS*B1];
    M3=[FS*A3,FS*B3];
    M=[M1;M3];
    
    end
    
    
    %%%%%%%%%% visualize CG result
    if visualization
    figure(21);
    subplot(311);     imshow(rr2(abs(first_guess)),[]); title('pinv original spectrum')
    subplot(312);     imshow(rr2(abs(pinvsol)),[]); title('pinv optimized spectrum')
    subplot(313);     imshow(rr2(abs(RCG)),[]); title('CG optimized spectrum')
    drawnow; colormap jet;
    end
    
    pause(1)
    
    end
    
end

fprintf('t: %i seconds per slice \n \n',toc)

end


%% TODO:  crop to original FOV, increase resolution, and normalize

ifft3_meas = @(I,dim) fftshift(ifftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
fft3_meas = @(I,dim) fftshift(fftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??

Recon_originalFOV=cell(2,2)
for ii=1:2;
    for jj=1:2
        if doubleFOV
            Recon_originalFOV{ii,jj}=Recon{ii,jj}(33:96,33:96,33:96);
            Recon_originalFOV{ii,jj}=fft3_meas(Recon_originalFOV{ii,jj});
            Recon_originalFOV{ii,jj}=padarray(Recon_originalFOV{ii,jj},[32 32 32]);
            Recon_originalFOV{ii,jj}=ifft3_meas(Recon_originalFOV{ii,jj});
            
        else
            Recon_originalFOV{ii,jj}=Recon{ii,jj};
            if ~zerofilloption
                Recon_originalFOV{ii,jj}=fft3_meas(Recon_originalFOV{ii,jj});
                Recon_originalFOV{ii,jj}=padarray(Recon_originalFOV{ii,jj},[32 32 32]);
                Recon_originalFOV{ii,jj}=ifft3_meas(Recon_originalFOV{ii,jj});
            end
        end
    end;
    
end

% imagine(Recon_originalFOV{1,1},Recon_originalFOV{1,2},Recon_originalFOV{2,1},Recon_originalFOV{2,2},protonimage)
% imagine(Recon_originalFOV{2,1},Recon_originalFOV{2,2},protonimage)
% imagine(Recon_originalFOV{1,1},Recon_originalFOV{1,2},protonimage)

%% OverLay Image 

params = CreateOverlayImageParams() 
% params.slicerange=65
params.n=1
params.export=0
params.fignumber=50
params.threshold=0.04
params.Alpha=0.3;

params.slicerange=70
params.maxfactor=20;
params.minfactor=3;
params.orientation='axi'
params.exportname='0h_PFOB_'
CreateOverlayImage(protonimage,squeeze(Recon_originalFOV{1,1}),[],params)

params.threshold=0.03
params.maxfactor=6;
params.minfactor=200;
params.fignumber=51
params.exportname='0h_PFCE_'
CreateOverlayImage(protonimage,[],squeeze(Recon_originalFOV{1,2}),params) 
%%
cd('L:\basic\divi\Projects\cosart\fluor\fluor\experiments\comparison_methods_01152018')

save('Comparison_005.mat','Recon_originalFOV')
