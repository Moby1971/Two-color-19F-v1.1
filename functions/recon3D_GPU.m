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

% folder='L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive\071118_19Fphantom_Oj1\';
folder='L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive 2\CC071318_19FPhantom2_Ol1'
protonimages='improton.mat'
kspaces='kspaces.mat'


doubleFOV=0;        % 1: double FOV in phase encoding dims/ 0: half FOV in readout dir
zerofilloption=0 

fprintf('F-19 reconstruction... \n')
cd(folder)
fprintf('loading kspace data... \n')
load(protonimages);
load(kspaces);

directions={'FH','HF','LR','RL'}

% preprocessing the data
run fix_double_FOV_4readouts.m

% initialize reconstruction cell
Recon{1,1}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{1,2}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{2,1}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{2,2}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{3,1}=zeros(size(K1,1),size(K1,2),size(K1,3));
Recon{3,2}=zeros(size(K1,1),size(K1,2),size(K1,3));

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

visualization=0; 

% parameters 
BW=4.4643e+04;
ppm=282.5685;
%%
fprintf('Define functions, operators and derive parameters from data...\n')

[nx,ny,nz]=size(K1); % imsizes
N=nx*ny*nz; % number of pixels in one image

if doubleFOV 
   BWpix=BW/(nx/2); % for a double FOV, the BW/pix is essentially doubled 
else
   BWpix=BW/(nx);
end

%functions
vec= @(I) reshape(I,[numel(I), 1]);
rr2 = @(I) reshape(I,[nx,ny,nz,2]);

rr1 = @(I) reshape(I,[nx,ny,nz,4]);
rr = @(I) reshape(I,[nx,ny,nz]);

% operatorsFSp
F3=opDFT3(nx,ny,nz,1);
FB=opBlockDiag(F3,F3,F3,F3);  % 4 fourier ops - one for each image 
TVOP2=MakeTVOperator(nx);

run spectrum_based_BW.m  %calculate PFOB spectrum 


%%
slice= 1; %slice to visualize 

data=[vec(K1); vec(K2); vec(K3); vec(K4)]; % data matrix 

if phaseremoval
fprintf('removing phase of image data... \n')
data= FB*abs(opInverse(FB)*data); % 
end

k1=data(1:N);
k2=data(N+1:2*N);
k3=data(2*N+1:3*N);
k4=data(3*N+1:4*N);

linear_recon=opInverse(FB)*data;    %linear recon of data 

%%%%%%%%%% visualize linear reconstruction 
if visualization; 
linear_recon_reshaped=rr1(linear_recon); linear_recon_reshaped=squeeze(linear_recon_reshaped(:,:,slice,1));
figure(1);subplot(211); imshow(abs(linear_recon_reshaped),[]); title('magnitude linear recon') 
figure(1);subplot(212); imshow(angle(linear_recon_reshaped),[-pi, pi]); title('phase linear recon');
end
imagine(rr1(linear_recon))
%% construct the measurement operator 
for ii=1:4; Spectrum{ii}=Spectrum1; end
M= MakeMeasurementOperator(Spectrum,nx,ny,nz,directions);
data=[k1;k2;k3;k4];

%% no l1-minimzation (first guess) 

first_guess=pinv(M)*data;
Recon{2,1}=rr(abs(first_guess(1:N)));
Recon{2,2}=rr(abs(first_guess(N+1:2*N)));

imagine(rr2(abs(first_guess)),'Name','Test')
%% l1-minimization CG 

lambda=25

if use_CG
    fprintf('CG recon... \n')

    for outeriter=1:3
        if outeriter==1; RCG=zeros(size(first_guess)); end
        RCG=nl_conjgrad_fluor_3D(M, data, RCG,niter,TVOP2,lambda,[nx ny nz],visualization,54);
    end
    Recon{1,1}=rr(abs(RCG(1:N)));
    Recon{1,2}=rr(abs(RCG(N+1:2*N)));
end

imagine(Recon)
%% Blind Deconvolution 
M= MakeMeasurementOperator(Spectrum,nx,ny,nz,directions);
if use_CG
    fprintf('Blind Deconvolution CG recon... \n')
    
    for outeriter=1:10
        fprintf('Outeriter %2d \n',outeriter)
        if outeriter==1; RCG=zeros(size(first_guess)); end
        RCG=nl_conjgrad_fluor_3D(M, data, RCG,niter,TVOP2,lambda,[nx ny nz],visualization,54);
        k = BlindDeconvolutionOptimization(RCG,N,linear_recon,ny,Spectrum1,rr,directions,1,visualization);
        
        UpdatedSpectrum=cell(1,4);
        UpdatedSpectrum{1}=k(1:nx).';
        UpdatedSpectrum{2}=k(nx+1:2*nx).';%k(nx+1:2*nx);
        UpdatedSpectrum{3}=k(2*nx+1:3*nx).';%k(2*nx+1:3*nx);
        UpdatedSpectrum{4}=k(3*nx+1:4*nx).';%k(3*nx+1:4*nx);

        M= MakeMeasurementOperator(UpdatedSpectrum,nx,ny,nz,directions);


    end
    Recon{3,1}=rr(abs(RCG(1:N)));
    Recon{3,2}=rr(abs(RCG(N+1:2*N)));
end

imagine(Recon)

%%
Recon{1,1}=fftshift(Recon{1,1})
Recon{1,2}=fftshift(Recon{1,2})
Recon{2,1}=fftshift(Recon{2,1})
Recon{2,2}=fftshift(Recon{2,2})
Recon{3,1}=fftshift(Recon{3,1})
Recon{3,2}=fftshift(Recon{3,2})

%%
ifft3_meas = @(I,dim) fftshift(ifftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
fft3_meas = @(I,dim) fftshift(fftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??

Recon_originalFOV=cell(3,2)
for ii=1:3;
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
save('Recon_scan2_17-7.mat','Recon_originalFOV')


%%

reconset=3

exportfolder='L:\basic\divi\Projects\cosart\fluor\fluor\recon_realdata\jul17\Recons\'
mkdir(exportfolder); cd(exportfolder)
params.orientation='cor'
params.fignumber=53
params.export=1
params.slicerange=45;
params.exportname='july18_2_1'


params.maxfactor=1
params.threshold=0.07
params.minfactor=1
params.n=1
CreateOverlayImage(improton,squeeze(Recon_originalFOV{reconset,1}),[],params)

params.exportname='july18_2_2'

params.minfactor=1;
params.maxfactor=1.8;
params.fignumber=54;

CreateOverlayImage(improton,[],squeeze(Recon_originalFOV{reconset,2}),params)




