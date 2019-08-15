function I= extract_slice (I,sl,dim)

if dim==3;
    I=squeeze(I(:,:,sl));
elseif dim==2;
    I=squeeze(I(:,:,sl));
elseif dim==1
    I=squeeze(I(sl,:,:));
else
    error('dim not ok')
end;
