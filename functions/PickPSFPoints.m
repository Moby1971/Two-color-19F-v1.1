function Spectrum=PickPSFPoints(l1,varargin)
nx=size(l1,1); 
BW=4.4643e+04 
BWpix=BW/nx;
ppm=282.5685;

if nargin==3
   saveoption=varargin{1} ;
   savename=varargin{2};
else; saveoption=0; end

figure(500);clf; imshow(abs(l1),[]);title('pick left upper and right lower points')
[xROI,yROI]=ginput(2);

l1zoom=l1(yROI(1):yROI(2),xROI(1):xROI(2));

figure(500);clf; 

imshow(abs(l1zoom),[]);title('pick the 5 PSF points (good order)')

[xPFOB,yPFOB]=ginput(5);

xPFOB=xPFOB+xROI(1);
yPFOB=yPFOB+yROI(1);

for ii=1:5
PFOBalpha1(ii)=l1(round(yPFOB(ii)),round(xPFOB(ii)));
end


[~,~,PFOB,~]=calcspectra_BW(ppm,BWpix); %126.4??
pixlocs=1+round(-PFOB); 
pixlocs(pixlocs<0)=nx+pixlocs(pixlocs<0);

PFOB_alpha=PFOBalpha1.';
Spectrum=zeros(1,nx);
Spectrum(pixlocs)=PFOB_alpha./sum(abs(PFOB_alpha(:)));
 if saveoption==1
     save([savename,'.mat'],'Spectrum'); 
 end

end