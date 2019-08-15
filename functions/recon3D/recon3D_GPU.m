% new attempt at 3D recon of fluor data 
% with GPU implementation 
% J SCHOORMANS 
% 5 january 2018 
clear all; close all; clc;
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\imagine'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor'))


%%

folder='L:\basic\divi\Projects\cosart\fluor\fluor\scans\jul17';
protonimages='protonim.mat'
kspaces='kspaces.mat'


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

%%
% reconstruction options:
dim= 3;             % shared phase_encoding dimension(??) 
phaseremoval=1;     % option to remove phase 
phasecheckerboardremoval=~phaseremoval; 

use_CG=1;           %use CG recon
lambda=70;
niter=15;

senselikePFCE=0;
normthreshold=.99;

use_empirical_alpha=1; 
% empirical_alpha=[64 101 64 175 64]; % peak heights PFOB
% empirical_alpha=[184 502 235 427 274];
empirical_alpha=[1.2 1.77 0.917 2.165 0.917]

visualization=1; 

% parameters 
BW=4.4643e+04;
ppm=282.5685;
%%
fprintf('Define functions, operators and derive parameters from data...\n')

nx=size(K1,1); % imsize 1
ny=size(K1,2); % imsize 2
nz=size(K1,3); % imsize 2

N=nx*ny*nz; % number of pixels in one image

if doubleFOV 
   BWpix=BW/(nx/2); % for a double FOV, the BW/pix is essentially doubled 
else
   BWpix=BW/(nx);
end

%functions
vec= @(I) reshape(I,[numel(I), 1]);
rr1 = @(I) reshape(I,[nx,ny,nz,2]);
rr = @(I) reshape(I,[nx,ny,nz]);

% operatorsFSp
F3=opDFT3(nx,ny,nz,1);
FB=opBlockDiag(F3,F3);  % 4 fourier ops - one for each image 

run spectrum_based_BW.m  %calculate PFOB spectrum 


%%
slice= 1; %slice to visualize 

data=[vec(K1); vec(K3)]; % data matrix 

if phaseremoval
fprintf('removing phase of image data... \n')
data= FB*abs(opInverse(FB)*data); % 
end

% if phasecheckerboardremoval
% fprintf('removing phase checkerboard image data... \n')
% data= (opInverse(F3)*data); % 
% checkerboard1=(-1).^linspace(1,ny,ny).*1i;
% checkerboard2=(-1).^linspace(1,nx*2,nx*2).*1i;
% data=rr(data);
% data=bsxfun(@times,data,checkerboard1.');
% data=bsxfun(@times,data,checkerboard2);
% data=F3*data(:);
% end

k1=data(1:N);
k3=data(N+1:2*N);

linear_recon=opInverse(FB)*data;    %linear recon of data 

%%%%%%%%%% visualize linear reconstruction 
if visualization; 
linear_recon_reshaped=rr1(linear_recon); linear_recon_reshaped=squeeze(linear_recon_reshaped(:,:,slice,1));
figure(1);subplot(211); imshow(abs(linear_recon_reshaped),[]); title('magnitude linear recon') 
figure(1);subplot(212); imshow(angle(linear_recon_reshaped),[-pi, pi]); title('phase linear recon');
end

% % reshaped magnitude of linear images
% l1=(rr1(linear_recon(1:N)));
% l3=(rr1(linear_recon((N)+1:2*(N))));


% make convolution operators (include fftshift in operators for simplicity)
Spectrum_vert=Spectrum1.';
Spectrum=(Spectrum1);
B=opDirac(nx*ny*nz); 
A1=opConvolve3D(nx,ny,nz,Spectrum.',[1 1 1],'cyclic');
A3=opConvolve3D(nx,ny,nz,Spectrum,[1 1 1],'cyclic');    %perhaps 65?


if senselikePFCE
%     fprintf('Recon: Sense-like PFCE (MORE TESTING NEEDED!) /n')
%     
%     l1norm=abs(l1(:));
%     l3norm=abs(l3(:));
%     sos=sqrt(l1norm.^2+l3norm.^2+1e-6);
%     allnorms=cat(3,l1norm./sos,l3norm./sos);
%     normstochange=(max(allnorms,[],3)>normthreshold);
%     
%     if visualization
%     figure(22);imshow(rr1(normstochange)); end
% 
%     allnorms(normstochange,:,:)=ones(sum(normstochange),1,2)./2;
%     sos=sqrt(sum((allnorms.^2),3)+1e-6);%recalculate SoS
%     
%     % to do : find out which of these options is the best
%     % B1=opDiag(allnorms(:,1,1)./sos);
%     % B3=opDiag(allnorms(:,1,2)./sos);
%     B1=opDiag(2*l1./(l1+l3));
%     B3=opDiag(2*l3./(l1+l3));
else
    B1=B;  B3=B;
end

% image-space recon ('first-guess')
M=[F3*A1,F3*B;F3*A3,F3*B];
first_guess=pinv(M)*data;
Recon{2,1}=rr(abs(first_guess(1:N)));
Recon{2,2}=rr(abs(first_guess(N+1:2*N)));


if use_CG
    fprintf('CG recon... \n')
    
    TVconvolution=zeros(nx,1,1);
    TVconvolution(1)=1; TVconvolution(end)=-1;
    TVx=opConvolve3D(nx,nx,nx,TVconvolution,[1 1 1],'cyclic');
    TVy=opConvolve3D(nx,nx,nx,permute(TVconvolution,[2 1 3]),[1 1 1],'cyclic');
    TVz=opConvolve3D(nx,nx,nx,permute(TVconvolution,[2 3 1]),[1 1 1],'cyclic');
    TVOP=[TVx;TVy;TVz];
    TVOP2=opBlockDiag(TVOP,TVOP)
    
    for outeriter=1:3
        if outeriter==1; RCG=zeros(size(first_guess)); end
        RCG=nl_conjgrad_fluor_3D(M, data, RCG,niter,TVOP2,lambda,[nx ny nz],1,3);
    end
    Recon{1,1}=rr(abs(RCG(1:N)));
    Recon{1,2}=rr(abs(RCG(N+1:2*N)));
end


if visualization
figure(5000); 
subplot(221); imshow(abs(squeeze(Recon{1,1}(:,:,slice))),[])
subplot(222); imshow(abs(squeeze(Recon{1,2}(:,:,slice))),[])
subplot(223); imshow(abs(squeeze(Recon{2,1}(:,:,slice))),[])
subplot(224); imshow(abs(squeeze(Recon{2,2}(:,:,slice))),[]); colormap jet
end

%%
Recon{1,1}=fftshift(Recon{1,1})
Recon{1,2}=fftshift(Recon{1,2})
%%
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

%%

exportfolder='L:\basic\divi\Projects\cosart\fluor\fluor\recon_realdata\dec17\Recons\Mouse4_10_1_2018'
mkdir(exportfolder); cd(exportfolder)
params.orientation='axi'
params.fignumber=53
params.export=1
params.slicerange=85;
params.exportname='Mouse2_48H_PFOB_empPSF'


params.maxfactor=19;
params.threshold=0.05
params.minfactor=1
params.n=1
CreateOverlayImage(protonimage,squeeze(Recon_originalFOV{1,1}),[],params)

params.exportname='Mouse2_48H_PFCE_empPSF'

params.minfactor=1;
params.maxfactor=1.8;
params.fignumber=54;

CreateOverlayImage(protonimage,[],squeeze(Recon_originalFOV{1,2}),params)




