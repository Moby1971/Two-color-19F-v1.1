function K1 = registration_correction_PFCE(app,K1,K3,visualization)
% outputs corrected K1 (translated) versus K3

L1=abs(ifftn(K1));
L3=abs(ifftn(K3));

[optimizer, metric] = imregconfig('monomodal');

[R_reg]=imregtform(L1,L3,'translation',optimizer, metric,'DisplayOptimization',0);
[dummy]=imregister(L1,L3,'translation',optimizer, metric,'DisplayOptimization',0);

shifts=R_reg.T(4,1:3);

TextMessage(app,strcat("Shifts: x = ",num2str(shifts(1)),"   y = ",num2str(shifts(2)),"   z = ",num2str(shifts(3))));

% do the correction in k-space (on L1?)
corrx=exp(-linspace(-pi,pi,size(K1,1)).*1i.*shifts(2));
corrx=permute(corrx,[2 1]);
corry=exp(-linspace(-pi,pi,size(K1,2)).*1i.*shifts(1));
corrz=exp(-linspace(-pi,pi,size(K1,3)).*1i.*shifts(3));
corrz=permute(corrz,[1 3 2]);
K1_corr=K1;
K1_corr=bsxfun(@times,K1_corr,corrx);
K1_corr=bsxfun(@times,K1_corr,corry);
K1_corr=bsxfun(@times,K1_corr,corrz);

% check correction
% L1_corr=abs(ifftn(K1_corr));
%if visualization
%imagine(L1_corr,'Name','L1corr',L1,'Name','L1 original',dummy,'Name','registered L1',L3,'Name','L3'); end
%apply correction

K1=K1_corr;