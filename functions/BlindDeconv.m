function Spectrum=BlindDeconv(l1,varargin)
nx=size(l1,1); 
BW=4.4643e+04 
BWpix=BW/nx;
ppm=282.5685;

if nargin==3
   saveoption=varargin{1} ;
   savename=varargin{2};
else; saveoption=0; end

figure(500);clf; imshow(abs(l1),[]);title('pick ROI')
[xROI,yROI]=ginput(2);
l1zoom=l1(yROI(1):yROI(2),xROI(1):xROI(2));


if size(l1zoom,1)>size(l1zoom,2)
    l1zoom=l1(:,xROI(1):xROI(2));
else
    l1zoom=l1(yROI(1):yROI(2),:);
end
figure(500);clf; 
imshow(abs(l1zoom),[]);title('Blind Deconv Region')

% calculate theoretical PSF
[~,~,PFOB,PFOB_alpha]=calcspectra_BW(ppm,BWpix); %126.4??
pixlocs=1+round(-PFOB); 
pixlocs(pixlocs<0)=nx+pixlocs(pixlocs<0);
Spectrum=zeros(1,nx);
Spectrum(pixlocs)=PFOB_alpha./sum(abs(PFOB_alpha(:)));
deconvreal = deconvreg(real(l1zoom),flip(Spectrum.',2),0);
deconvimag = deconvreg(imag(l1zoom),Spectrum.',0);


figure(501); subplot(1,2,1); imshow(abs(cat(2,real(l1zoom),deconvreal)),[])
figure(501); subplot(1,2,2); imshow(abs(cat(2,imag(l1zoom),deconvimag)),[])



 if saveoption==1
     save([savename,'.mat'],'Spectrum'); 
 end

end