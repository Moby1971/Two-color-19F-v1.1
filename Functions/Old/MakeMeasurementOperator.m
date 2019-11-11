function M= MakeMeasurementOperator(Spectrum,nx,ny,nz,directions,varargin)
% Constructs the Matrix (F 0/0 F)(C 0/0 C)

nacq=length(directions); % number of acquisitions in system
assert(length(Spectrum)==nacq,'Spectrum should be a cell'); %Spectrum should be a cell

if nargin>5
    F=varargin{1};
else
    for ii=1:nacq
        F{ii}=opDFT3(nx,ny,nz,1); end
end

for i=1:nacq 
assert(isrow(Spectrum{1}),'Spectrum should be a row vector')
end

B=MakeIdentityOperator(nx,ny,nz);
for ii=1:nacq
    S=MakeSpectrum(Spectrum{ii},directions{ii});
    A{ii}=MakeConvOperator(nx,ny,nz,S);
end

M=[];
for ii=1:nacq
    M=[M;F{ii}*A{ii},F{ii}*B];
end
end

function A=MakeConvOperator(nx,ny,nz,S)
A=opConvolve3D(nx,ny,nz,S,[1 1 1],'cyclic');
end

function B=MakeIdentityOperator(nx,ny,nz)
B=opDirac(nx*ny*nz);
end

function S=MakeSpectrum(Spectrum,dir)
if dir=='HF'
    S=Spectrum;
elseif dir=='FH'
    S=flip(Spectrum);
elseif dir=='LR'
    S=Spectrum.';
elseif dir=='RL'
    S=flip(Spectrum.');
else
    error('direction error')
end
end
