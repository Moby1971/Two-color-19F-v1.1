function TVOP2=MakeTVOperator(nx)
% Total Variation Operator (Three-dimensional)

TVconvolution=zeros(nx,1,1);
TVconvolution(1)=1; TVconvolution(end)=-1;
TVx=opConvolve3D(nx,nx,nx,TVconvolution,[1 1 1],'cyclic');
TVy=opConvolve3D(nx,nx,nx,permute(TVconvolution,[2 1 3]),[1 1 1],'cyclic');
TVz=opConvolve3D(nx,nx,nx,permute(TVconvolution,[2 3 1]),[1 1 1],'cyclic');
TVOP=[TVx;TVy;TVz];
TVOP2=opBlockDiag(TVOP,TVOP);

end