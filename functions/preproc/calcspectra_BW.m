function [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(ppm_hz,BW_pix)
%inputs:
% ppm_hz : ppm in Hz
% BW_pix : BW per pixel

%outputs: 
% PFCE,PFCE_alpha,PFOB,PFOB_alpha:locations (pixels) and intensities (A.U.) of chemical shift peaks  


%% FOR PFCE
PFCE=0;
PFCE_alpha=20;

%% FOR PFOB
df=[26.7,8.6,-27.1,-31.6,-36.2]; %ppm

dx=(ppm_hz.*df)./BW_pix;
PFOB=dx; % in pixels
PFOB_alpha=[2,3,2,6,2];

end