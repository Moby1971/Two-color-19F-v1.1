% test 3D operators SPOT 
% J SCHOORMANS 9 JNAUARY 2018



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
%% options 
nx=64
ny=64
nz=1
useGPU=0
ignoreFFT=1

rr= @(x) reshape(x,[nx nx nz 2]); %reshape function 
rr2= @(x) x(:); % reshape function 
N=(nx*ny*nz*2); 
Data=cat(4,K1,K3);  %4d matrix (kx,ky,kz,direction)
Data=Data(:); 

I_input=zeros(nx,nx,1,2); 
I_input(49,30:38,:,1)=ones(1,9,size(I_input,3),1); 
I_input(20:23,20:23,:,2)=ones(4,4,size(I_input,3),1); 

%% SPOT FUNCTION
%%%% operator from image space (PFCE/PFOB) to measured image space (direction 1/direction 2) 
%[I, C_H; I, C_V]

% three dimensional spectra
SpectrumPFOB1=zeros(64,64,nz); SpectrumPFOB1(:,1,1)=Spectrum; 
SpectrumPFOB2=zeros(64,64,nz); SpectrumPFOB2(1,:,1)=Spectrum;  

Conv_op= @(x,b) ifftn(fftn(x).*fftn(b)); 
Forward_op= @(x) cat(4,x(:,:,:,2),x(:,:,:,2)) +...
    cat(4,Conv_op(x(:,:,:,1),SpectrumPFOB1),Conv_op(x(:,:,:,1),SpectrumPFOB2))

Measurement_Forward = @(x,mode) rr2((Forward_op(rr(x)))); 
A=opFunction(N,N, Measurement_Forward);

%%==========================================

I_measured=Measurement_Forward((I_input));
I_measured=(rr(I_measured));
Data_2=rr2((I_measured)); 
first_guess=(rr(A\Data_2));

figure(10); 
subplot(321); imshow(abs(squeeze(I_input(:,:,1,1))),[]); title('Data 1 (PFOB)')
subplot(322); imshow(abs(squeeze(I_input(:,:,1,2))),[]); title('Data 1 (PFCE)')
subplot(323); imshow(abs(squeeze(I_measured(:,:,1,1))),[]); title('A*Data')
subplot(324); imshow(abs(squeeze(I_measured(:,:,1,2))),[]); title('A*Data')
subplot(325); imshow(abs(squeeze(first_guess(:,:,nz,1))),[]); title('A\\(A*Data)')
subplot(326); imshow(abs(squeeze(first_guess(:,:,nz,2))),[]); title('A\\(A*Data)'); colormap jet

%% SPOT FUNCTION 2D (SAME PROBLEMS)

rr= @(x) reshape(x,[nx ny 2]); %reshape function 
rr2= @(x) x(:); % reshape function 

SpectrumPFOB1=zeros(64,64); SpectrumPFOB1(:,1)=Spectrum; 
SpectrumPFOB2=zeros(64,64); SpectrumPFOB2(1,:)=Spectrum;  

Conv_op= @(x,b) ifftn(fftn(x).*fftn(b)); 
Forward_op= @(x) cat(3,x(:,:,2),x(:,:,2)) +...
    cat(3,Conv_op(x(:,:,1),SpectrumPFOB1),Conv_op(x(:,:,1),SpectrumPFOB2))

Measurement_Forward = @(x,mode) rr2((Forward_op(rr(x)))); 
A=opFunction(N,N, Measurement_Forward);

%%==========================================

I_measured=Measurement_Forward((I_input));
I_measured=(rr(I_measured));
Data_2=rr2((I_measured)); 
first_guess=(rr(A\Data_2));

figure(12); 
subplot(321); imshow(abs(squeeze(I_input(:,:,1))),[]); title('Data 1 (PFOB)')
subplot(322); imshow(abs(squeeze(I_input(:,:,2))),[]); title('Data 1 (PFCE)')
subplot(323); imshow(abs(squeeze(I_measured(:,:,1))),[]); title('A*Data')
subplot(324); imshow(abs(squeeze(I_measured(:,:,2))),[]); title('A*Data')
subplot(325); imshow(abs(squeeze(first_guess(:,:,1))),[]); title('A\\(A*Data)')
subplot(326); imshow(abs(squeeze(first_guess(:,:,2))),[]); title('A\\(A*Data)'); colormap jet



%% SPOT 2D OPERATORS (NEGLECTING FFT FOR NOW!)
A1=opConvolve(nx,ny,Spectrum.',[1 1],'cyclic');
A3=opConvolve(nx,ny,Spectrum,[1 1],'cyclic'); %perhaps 65?
B=opDirac(nx*ny); 

M=[A1,B;A3,B]

y_SPOT=rr(M*I_input(:));

x_hat_SPOT=rr(M\y_SPOT(:));

figure(11); 
subplot(321); imshow(abs(squeeze(I_input(:,:,1,1))),[]); title('Data 1  (PFOB)')
subplot(322); imshow(abs(squeeze(I_input(:,:,1,2))),[]); title('Data 1  (PFCE)')
subplot(323); imshow(abs(squeeze(y_SPOT(:,:,1,1))),[]); title('M*Data')
subplot(324); imshow(abs(squeeze(y_SPOT(:,:,1,2))),[]); title('M*Data')
subplot(325); imshow(abs(squeeze(x_hat_SPOT(:,:,nz,1))),[]); title('A\\(A*Data)')
subplot(326); imshow(abs(squeeze(x_hat_SPOT(:,:,nz,2))),[]); title('A\\(A*Data)') ;colormap jet;

%% 
%=================================================================
%=================================================================
%=================================================================
%=================================================================
%=============   NEW 3D CONVOLUTION OPERATOR  ====================
%=================================================================
%=================================================================
%=================================================================
%=================================================================

A=opConvolve3D(64,64,64,rand(1,1,64),[1,1,1],'cyclic')

x=zeros(64,64,64); x(32,32,32)=1; 
y=A*x(:);
y=reshape(y,[64,64,64]); 

figure(1); 
subplot(131); imshow(abs(squeeze(y(:,:,32))),[]);
subplot(132); imshow(abs(squeeze(y(:,32,:))),[]);
subplot(133); imshow(abs(squeeze(y(32,:,:))),[]);
colormap jet

xhat=A\y(:); 
xhat=reshape(xhat,[64 64 64]); 

figure(2); 
subplot(131); imshow(abs(squeeze(xhat(:,:,32))),[]);
subplot(132); imshow(abs(squeeze(xhat(:,32,:))),[]);
subplot(133); imshow(abs(squeeze(xhat(32,:,:))),[]);
colormap jet

%% testing of newly made opConvolve3D operator 
nz=64; ny=64; nz=64; 

I_input=zeros(nx,nx,nz,2); 
I_input(49,30:38,:,1)=ones(1,9,size(I_input,3),1); 
I_input(20:23,20:23,:,2)=ones(4,4,size(I_input,3),1); 
rr= @(x) reshape(x,[nx nx nz 2]); %reshape function 

A1=opConvolve3D(nx,ny,nz,Spectrum.',[1 1 1],'cyclic');
A3=opConvolve3D(nx,ny,nz,Spectrum,[1 1 1],'cyclic'); %perhaps 65?
B=opDirac(nx*ny*nz); 
M=[A1,B;A3,B]

y_SPOT=rr(M*I_input(:));
x_hat_SPOT=rr(M\y_SPOT(:));

figure(15); 
subplot(321); imshow(abs(squeeze(I_input(:,:,1,1))),[]); title('Data 1  (PFOB)')
subplot(322); imshow(abs(squeeze(I_input(:,:,1,2))),[]); title('Data 1  (PFCE)')
subplot(323); imshow(abs(squeeze(y_SPOT(:,:,1,1))),[]); title('M*Data')
subplot(324); imshow(abs(squeeze(y_SPOT(:,:,1,2))),[]); title('M*Data')
subplot(325); imshow(abs(squeeze(x_hat_SPOT(:,:,nz,1))),[]); title('A\\(A*Data)')
subplot(326); imshow(abs(squeeze(x_hat_SPOT(:,:,nz,2))),[]); title('A\\(A*Data)') ;colormap jet;

%%

%=================================================================
%=================================================================
%=================================================================
%=================================================================
%=================   NEW 3D FFT OPERATOR  ========================
%=================================================================
%=================================================================
%=================================================================
%=================================================================

F3=opDFT3(64,64,64,1)
x=repmat(phantom(64),[1 1 64])
y=F3*x(:);
xhat=F3'*y;
plot(xhat,x(:),'');

xhat=reshape(xhat,[64,64,64]);
figure(16); 
subplot(231); imshow(abs(squeeze(x(:,:,32))),[]);
subplot(232); imshow(abs(squeeze(x(:,32,:))),[]);
subplot(233); imshow(abs(squeeze(x(32,:,:))),[]);
subplot(234); imshow(abs(squeeze(xhat(:,:,1))),[]);
subplot(235); imshow(abs(squeeze(xhat(:,32,:))),[]);
subplot(236); imshow(abs(squeeze(xhat(32,:,:))),[]);
colormap jet

%% Combine New 3D operators!!!

M=[F3*A1,F3*B;F3*A3,F3*B]

y_SPOT=rr(M*I_input(:));
x_hat_SPOT=rr(M\y_SPOT(:));


figure(17); 
subplot(321); imshow(abs(squeeze(I_input(:,:,32,1))),[]); title('Data 1  (PFOB)')
subplot(322); imshow(abs(squeeze(I_input(:,:,32,2))),[]); title('Data 1  (PFCE)')
subplot(323); imshow(abs(squeeze(y_SPOT(:,:,nz,1))),[]); title('M*Data')
subplot(324); imshow(abs(squeeze(y_SPOT(:,:,nz,2))),[]); title('M*Data')
subplot(325); imshow(abs(squeeze(x_hat_SPOT(:,:,32,1))),[]); title('A\\(A*Data)')
subplot(326); imshow(abs(squeeze(x_hat_SPOT(:,:,32,2))),[]); title('A\\(A*Data)') ;colormap jet;

%%
%=================================================================
%=================================================================
%=================================================================
%=================================================================
%==========  USE 3D FFT OPERATOR WITH REAL DATA  =================
%=================================================================
%=================================================================
%=================================================================
%=================================================================

rr= @(x) reshape(x,[64 64 64 2]); %reshape function 



A1=opConvolve3D(64,64,64,Spectrum.',[1 1 1],'cyclic');
A3=opConvolve3D(64,64,64,Spectrum,[1 1 1],'cyclic');    %perhaps 65?
B=opDirac(64^3); 
F3=opDFT3(64,64,64,0);

M=[F3*A1,F3*B;F3*A3,F3*B];
Data=cat(4,K1,K3); Data=Data(:); 

F3_2=opBlockDiag(F3,F3);

Data= F3_2*abs(opInverse(F3_2)*Data); % 

Recon=rr(M\Data); 

RiFFT=rr(F3_2'*Data); 
%%
i=10

figure(17);
subplot(421); imshow(abs(squeeze(RiFFT(:,:,i,1))),[0 max(abs(RiFFT(:)))]); title('zerofill1')
subplot(422); imshow(abs(squeeze(RiFFT(:,:,i,2))),[0 max(abs(RiFFT(:)))]); title('zerofill2') ;colormap jet;
subplot(423); imshow(abs(squeeze(Recon(:,:,i,1))),[0 max(abs(Recon(:)))]); title('pinv 1')
subplot(424); imshow(abs(squeeze(Recon(:,:,i,2))),[0 max(abs(Recon(:)))]); title('pinv2') ;colormap jet;

%% 2D case 
dim=3
ifft_meas = @(I,dim) fftshift(ifft(ifftshift(I,dim),[],dim),dim); % iFFT in measurement direction  %%? which dimension??
vec= @(I) reshape(I,[numel(I), 1]);

sl=41

k1=extract_slice(ifft_meas(K1,dim),sl,3);
k3=extract_slice(ifft_meas(K3,dim),sl,3);
data=[vec(k1); vec(k3)]; % data matrix 

Spectrum_vert=Spectrum1.';
Spectrum=(Spectrum1);
A1=opConvolve(nx,ny,Spectrum_vert,[1 1],'cyclic');
A3=opConvolve(nx,ny,Spectrum1,[1 1],'cyclic'); %perhaps 65?
B=opDirac(nx*ny); 
ShiftOp=opConvolve(nx,ny,1,[nx/2 ny/2],'cyclic');
F=opDFT2(nx,ny,1); %fourier operator
FS=F*ShiftOp;

F2=opBlockDiag(FS,FS); 
data= F2*abs(opInverse(F2)*data); % 


M1=[FS*A1,FS*B];
M3=[FS*A3,FS*B];
M=[M1;M3];

first_guess=pinv(M)*data; first_guess=reshape(first_guess,[64 64 2]);
zf2d=reshape(F2'*data,[64 64 2]); 

subplot(425); imshow(abs(squeeze(zf2d(:,:,1))),[]); title('pinv 1'); title('2D Recon')
subplot(426); imshow(abs(squeeze(zf2d(:,:,1))),[]); title('pinv 1'); title('2D Recon')

subplot(427); imshow(abs(squeeze(first_guess(:,:,1))),[0 max(abs(first_guess(:)))]); title('pinv 1'); title('2D Recon')
subplot(428); imshow(abs(squeeze(first_guess(:,:,2))),[0 max(abs(first_guess(:)))]); title('pinv2') ;colormap jet;


%%



