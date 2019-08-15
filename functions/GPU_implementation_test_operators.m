% kladblok GPU implementatie FLUOR code 
clear all
clear all; close all; clc;
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\imagine'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor'))


folder='L:\basic\divi\Projects\cosart\fluor\fluor\recon_realdata\nov17';
cd(folder)
kspaces='kspaces112917_1.mat'
load(kspaces)
K1=K1(1:2:end,:,:);
K3=K3(:,1:2:end,:);

% parameters 
BW=4.4643e+04;
ppm=282.5685;
nx=64; visualization=0;
BWpix=BW/64;
use_empirical_alpha=0;
run spectrum_based_BW.m  %calculate PFOB spectrum 
Spectrum=Spectrum1;

%% 3D FFT
% test fluor recon operaties - 3D op de GPU 
K1G=gpuArray(K1); 

tic 
fftshift(fftn(ifftshift(K1)));
tFFT_CPU=toc; 


tic 
fftshift(fftn(ifftshift(K1G)));
tFFT_GPU=toc; 

fprintf('FFT: \n')
fprintf('Time taken CPU: %f \n',tFFT_CPU)
fprintf('Time taken GPU: %f \n',tFFT_GPU)
fprintf('GPU is %f  times faster \n\n',tFFT_CPU/tFFT_GPU)

%%% CONVN 
B=rand(1,size(K1,2));
BG=gpuArray(B);

tic 
convn(K1,B);
t_convn_CPU=toc;

tic
convn(K1G,BG);
t_convn_GPU=toc;

fprintf('convn: \n')
fprintf('Time taken CPU: %f \n',t_convn_CPU)
fprintf('Time taken GPU: %f \n',t_convn_GPU)
fprintf('GPU is %f  times faster \n\n',t_convn_CPU/t_convn_GPU)

%%
b=cat(4,K1,K3); %% 05-01-2018 copied from lower in the code

%%
iFFT_c_dim = @(x,dim) fftshift(ifft(ifftshift(x,dim),[],dim),dim);

SpectrumK1=(Spectrum.');
SpectrumK3=(Spectrum);
x0=iFFT_c_dim(iFFT_c_dim(iFFT_c_dim(b,1),2),3);
sl=24

figure(2); imshow(cat(2,abs(squeeze(x0(:,:,sl,1))),abs(squeeze(x0(:,:,sl,2)))),[]); drawnow; 
title('zerofill');
figure(3); imshow(cat(2,bsxfun(@times,ones(64,64),SpectrumK1),bsxfun(@times,ones(64,64),SpectrumK3)),[])
title('PSF with direction for meas 1 and meas 2')
%% 3D circular convolution in k- space - multiplication of kspaces; 
Conv_op_CPU= @(x,b) ifftn(fftn(x).*fftn(b)); 

%image space 3D PSF
SK1=zeros(64,64,64); SK1(1,:,1)=Spectrum.'; 
SK3=zeros(64,64,64); SK3(:,1,1)=Spectrum;  

x0=zeros(64,64,64); x0(1,1,1)=1;

x0_3dconvolved_1=(Conv_op_CPU(x0,SK1));
x0_3dconvolved_3=(Conv_op_CPU(x0,SK3));

figure(3); imshow(cat(2,squeeze(x0_3dconvolved_1(:,:,1)),squeeze(x0_3dconvolved_3(:,:,1))),[0 0.2]); 
title('circ 3d convolution -CPU')

%% 3D circular deconvolution 
Deconv_op_CPU= @(x,b) ifftn(fftn(x)./fftn(b)); 

x0_3d_deconvolved_1=(Deconv_op_CPU(x0_3dconvolved_1,SK1));
x0_3d_deconvolved_3=(Deconv_op_CPU(x0_3dconvolved_3,SK3));

figure(4);
imshow(cat(2,squeeze(x0_3d_deconvolved_1(:,:,1)),squeeze(x0_3d_deconvolved_3(:,:,1))),[0 0.2]);
title('circ 3d deconvolution - CPU')

%% compare circular convolution operator - CPU and GPU

B=rand([64,64,64]);
BG=gpuArray(B);
SK1G=gpuArray(SK1); 

tic ;
for i=1:500
temp=Conv_op_CPU(B,SK1);
temp=Conv_op_CPU(temp,SK1);
end
tCPU=toc;

tic ;
for i=1:500
temp=Conv_op_CPU(BG,SK1G);
temp=Conv_op_CPU(temp,SK1G);
end
tGPU=toc;

fprintf('circular convolution (500 times back/forth): \n')
fprintf('Time taken CPU: %f \n',tCPU)
fprintf('Time taken GPU: %f \n',tGPU)
fprintf('GPU is %f  times faster \n\n',tCPU/tGPU)




%% 3D CPU-operator based linear conjugate gradient algo

FFT_c_dim = @(x,dim) fftshift(fft(ifftshift(x,dim),[],dim),dim);
FFT_op_CPU= @(x) FFT_c_dim(FFT_c_dim(FFT_c_dim(x,1),2),3);


SpectrumFFT=FFT_op_CPU(cat(4,SK1,SK3));

%define forward operator
A= @(x) FFT_op_CPU(cat(4,x(:,:,:,2),x(:,:,:,2)))+FFT_op_CPU(cat(4,x(:,:,:,1),x(:,:,:,1))).*SpectrumFFT


%%%  start  algo 
niter=500;
b=cat(4,K1,K3);
x0=zeros(size(b));
x=x0;

r = b - A(x);
p = r;
rsold = r(:).' * r(:);
fprintf('rsold: %i \n', rsold)

for iter=1:niter
    Ap = A(p);
    alpha = rsold / (p(:).' * Ap(:));
    
    x = x + alpha * p;
    
    figure(1); imshow(cat(2,abs(squeeze(x(:,:,33,1))),abs(squeeze(x(:,:,33,2)))),[]); drawnow; 
    pause(1)
    
    r = r - alpha * Ap;
    rsnew = r(:).' * r(:);
    
    fprintf('iteration: %i      |rsnew: %i      |alpha: %f \n',iter,rsnew,alpha)
    if sqrt(rsnew) < 1e-10
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
end

%% if the fourier-convolution matrix (of 1d convolution) is a repeat of kx*kz the same k-line, 
% couldnt we find the PSF by an autocorrelation of the kspace 


for sl2=1:64
    sl2
    for sl=1:64
        kline=K1(:,sl,sl2);
        
        for shift=1:64;
            C(shift,sl,sl2)= corr(kline,circshift(kline,shift));
        end
    end
end
figure(5); plot(mean(mean(abs(C),2),3))


%%















