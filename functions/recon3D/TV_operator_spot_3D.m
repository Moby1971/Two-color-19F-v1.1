% total variation operator using 3D SPOT convolve 

TVconvolution=zeros(128,1,1); 
TVconvolution(1)=1; TVconvolution(end)=-1; 
TVx=opConvolve3D(128,128,128,TVconvolution,[1 1 1],'cyclic');
TVy=opConvolve3D(128,128,128,permute(TVconvolution,[2 1 3]),[1 1 1],'cyclic');
TVz=opConvolve3D(128,128,128,permute(TVconvolution,[2 3 1]),[1 1 1],'cyclic');
TVOP=[TVx;TVy;TVz];

%%
fprintf('\n\n\n')
size(phantomdata(:))

phantomTV=TVOP*phantomdata(:); 
size(phantomTV)

phantomhat=TVOP'*phantomTV(:);

size(phantomhat)
%%