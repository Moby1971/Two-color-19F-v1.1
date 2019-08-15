%% dealing with double FOV in readout direction (2 options)

if ~doubleFOV
    fprintf('Removing odd k-space points \n')
    K1=K1(1:2:end,:,:);
    K2=K2(1:2:end,:,:);
    K3=K3(:,1:2:end,:);
    K4=K4(:,1:2:end,:);
    
    if zerofilloption %double resolution
        K1=padarray(K1,[32 32 32]);
        K2=padarray(K2,[32 32 32]);
        K3=padarray(K3,[32 32 32]);
        K4=padarray(K4,[32 32 32]);
    end
    
    ifft3_meas = @(I,dim) fftshift(ifftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
    I1=ifft3_meas(K1);
    I2=ifft3_meas(K2);
    I3=ifft3_meas(K3);
    I4=ifft3_meas(K4);
    
    %   imagine(I1,I3)
    
end

if doubleFOV;
    fprintf('doubling FOV in two phase encoding dimensions... \n')
    
    fft3_meas = @(I,dim) fftshift(fftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
    ifft3_meas = @(I,dim) fftshift(ifftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
    
    I1=ifft3_meas(K1);
    I1=padarray(I1,[0 32 32]);
    K1=fft3_meas(I1);
    
    I3=ifft3_meas(K3);
    I3=padarray(I3,[32 0 32]);
    K3=fft3_meas(I3);
    
    I2=ifft3_meas(K2);
    I2=padarray(I2,[32 0 32]);
    K2=fft3_meas(I2);
    
    I4=ifft3_meas(K4);
    I4=padarray(I4,[32 0 32]);
    K4=fft3_meas(I4);
end
