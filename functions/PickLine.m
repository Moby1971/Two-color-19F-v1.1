
function Spectrum=PickLine(l1,mirror, varargin)
nx=size(l1,1); 
BW=4.4643e+04 
BWpix=BW/nx;
ppm=282.5685;

if nargin==4
   saveoption=varargin{1} ;
   savename=varargin{2};

elseif nargin==5
    coords=varargin{3}
else
    ; saveoption=0; 
end


figure(500);clf; imshow(abs(l1),[]);title('1: zoom pick line (begin/end)')
[xROI,yROI]=ginput(2);
l1zoom=l1(yROI(1):yROI(2),xROI(1):xROI(2));
figure(500);clf; imshow(abs(l1zoom),[]);
[xROI2,yROI2]=ginput(2);
xROI=xROI2+xROI(1);
yROI=yROI2+yROI(1); 

thr=@(Sp1) Sp1.*(abs(Sp1)>max(abs(Sp1(:)))/8);


if size(l1zoom,1)>size(l1zoom,2)
    linespectrum=l1(:,round(xROI(1)));
else
    linespectrum=l1(round(yROI(1)),:);
end



linespectrum=thr(linespectrum)
linespectrum=flip(linespectrum,1); 


% calculate theoretical PSF
[~,~,PFOB,PFOB_alpha]=calcspectra_BW(ppm,BWpix); %126.4??
pixlocs=1+round(-PFOB); 
pixlocs(pixlocs<0)=nx+pixlocs(pixlocs<0);
Spectrum=zeros(1,nx);
Spectrum(pixlocs)=PFOB_alpha./sum(abs(PFOB_alpha(:)));
SpectrumBW=Spectrum;
 if mirror; linespectrum=flip(linespectrum,1); end
linespectrum=linespectrum./sum(abs(linespectrum(:))); 
linespectrum=linespectrum.'; % to do: shift

figure(501);clf
hold on; 
plot(abs(linespectrum)); 
plot(abs(SpectrumBW)); 
hold off
[shiftx,shifty]=ginput(2)
shift=round(shiftx(1)-shiftx(2))
linespectrum=circshift(linespectrum,[0 shift])

figure(501);clf
hold on; 
plot(abs(linespectrum)); 
plot(abs(SpectrumBW)); 
hold off; pause(1)

Spectrum=linespectrum; 
if ~mirror 
    Spectrum=flip(Spectrum,1);
end
 if saveoption==1
     save([savename,'.mat'],'Spectrum'); 
 end

