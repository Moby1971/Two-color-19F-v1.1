function  CreateOverlayImage(protonimage,PFOB,PFCE,params) 

% 16 11 2017 J SCHOORMANS 
normalize = @(x) x./(median(abs(x(:)))+eps); 

if isempty(PFOB)&&isempty(PFCE);
    error('please input at least one overlay image')
else
    if isempty(PFOB);
        PFOB=zeros(size(PFCE));
    end
    if isempty(PFCE);
        PFCE=zeros(size(PFOB));
    end
end

% PFOB=normalize(PFOB); 
% PFCE=normalize(PFCE);
% protonimage=normalize(protonimage); 

%% params 
orientation=params.orientation ;
export=params.export; 
slicerange=params.slicerange;
thr = params.threshold; 
maxfactor=params.maxfactor ;
minfactor=params.minfactor;
exportname=params.exportname;
n=params.n; %scaling function
%% 

clear greenmap redmap
for ii=1:100; greenmap(:,ii)=[0 ii/100 0]; end; greenmap=greenmap.';
for ii=1:100; redmap(:,ii)=[ii/100 0 0]; end; redmap=redmap.';
colorthreshold = @(colorim,thr) colorim.*(colorim>thr); 
% colorthreshold = @(colorim,thr) colorim(colorim<thr)=0;  %all under threshold nan??

thrPFCE=abs(thr*max(PFCE(:)));
thrPFOB=abs(thr*max(PFOB(:)));

% maxcolorim=abs(max(PFOB(:)))/maxfactor;
mtemp=sort(abs(max(PFOB(:))))
maxcolorim=mtemp(round(0.9*length(mtemp)))/maxfactor;
mincolorim=abs(min(PFOB(:)))*minfactor;

% maxcolorim2=abs(max(PFCE(:)))/maxfactor;
mtemp=sort(abs(max(PFCE(:))))
maxcolorim2=mtemp(round(0.9*length(mtemp)))/maxfactor;
mincolorim2=abs(min(PFCE(:)))*minfactor;

maxprotonim=double(max(abs(protonimage(:))));

for sl=slicerange;
    switch orientation
        case 'cor'
            protonim=abs(squeeze(protonimage(:,:,sl)));
            protonim=mat2gray(protonim,[0 maxprotonim]);
            colorim=abs(squeeze(PFOB(:,:,sl)));
            colorim=colorthreshold(colorim,thrPFOB);
            colorim2=abs(squeeze(PFCE(:,:,sl)));
            colorim2=colorthreshold(colorim2,thrPFCE);
            
        case 'axi'
            protonim=abs(squeeze(protonimage(sl,:,:))).';
            protonim=mat2gray(protonim,[0 maxprotonim]);
            colorim=abs(squeeze(PFOB(sl,:,:))).';
            colorim=colorthreshold(colorim,thrPFOB);
            colorim2=abs(squeeze(PFCE(sl,:,:))).';
            colorim2=colorthreshold(colorim2,thrPFCE);
            
        case 'sag'
            protonim=abs(squeeze(protonimage(:,sl,:))).';
            protonim=mat2gray(protonim,[0 maxprotonim]);
            colorim=abs(squeeze(PFOB(:,sl,:))).';
            colorim=colorthreshold(colorim,thrPFOB);
            colorim2=abs(squeeze(PFCE(:,sl,:))).';
            colorim2=colorthreshold(colorim2,thrPFCE);
    end
 

fig3=figure(params.fignumber);clf;
axis off; set(fig3,'Color','White')

ax1 = axes('Parent',fig3);
ax2 = axes('Parent',fig3);
ax3 = axes('Parent',fig3);

set(ax1,'Visible','off');
set(ax2,'Visible','off');
set(ax3,'Visible','off');

if sum(colorim(:))>0
% h2=imshow(colorim,[mincolorim maxcolorim],'Parent',ax2);
% set(h2,'AlphaData',Alpha)

% new method makes [0,1] map --> scales with alphadata 
h2=imshow(colorim>0,[0 1],'Parent',ax2);
Alpha=(((abs(colorim)-abs(thrPFOB))./maxcolorim));
Alpha(Alpha<0)=0; Alpha(Alpha>1)=1; Alpha=Alpha.^n; 
set(h2,'AlphaData',Alpha)
colormap(ax2,greenmap)


end 

if sum(colorim2(:))>0
% h3=imshow(colorim2,[mincolorim2 maxcolorim2],'Parent',ax3);
% set(h3,'AlphaData',Alpha)

% new method makes [0,1] map --> scales with alphadata 
h3=imshow(colorim2>0,[0 1],'Parent',ax3);
Alpha=(((abs(colorim2)-abs(thrPFCE))./maxcolorim2));
Alpha(Alpha<0)=0; Alpha(Alpha>1)=1; Alpha=Alpha.^n; 
set(h3,'AlphaData',Alpha)
colormap(ax3,redmap)
end

h1=imshow(protonim,[0 1],'Parent',ax1);
colormap(ax1,'gray')

pause(0.1)
if export
eval(['export_fig ',exportname,orientation,'_',num2str(sl),'.tiff']); end
end

% figure(1000); subplot(212); plot(sort(Alpha(:))); 


end 



