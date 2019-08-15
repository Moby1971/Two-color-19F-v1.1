
fprintf('Calculating the spectrum of PFOB and PFCE \n')
%% calculating the spectrum of PFOB and PFCE
[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(ppm,BWpix); %126.4??
if use_empirical_alpha==true;
    PFOB_alpha=empirical_alpha; %based on l2 sl62, see notes
    if ~phaseremoval;
        PFOB_alpha=[4.65,6.57,-1.87,4.84,-3.34]+ 1i*[-1.61,-2.58,4.57,0.67,-3.79];
    end
end
pixlocs=1+round(-PFOB);
pixlocs(pixlocs<0)=nx+pixlocs(pixlocs<0);
Spectrum1=zeros(1,nx); Spectrum1(pixlocs)=PFOB_alpha./sum(PFOB_alpha(:));

%%
%%%%%%%%%% visualize spectrum
if visualization
    figure(101);
    plot(abs(Spectrum1)); hold on;
    title('Spectrum used in recon (?)')
end