function M = MakeMeasurementOperator2C(SpectrumC1,SpectrumC2,nx,ny,nz,directions,varargin)
% Constructs the Matrix (F 0/0 F)(C 0/0 C)

nacq=length(directions); % number of acquisitions in system

if nargin>6
    F = varargin{1};
else
    for ii = 1:nacq
        F{ii} = opDFT3(nx,ny,nz,1);
    end
end


% Compound 1
for ii=1:nacq
    S1 = MakeSpectrum(SpectrumC1{ii},directions{ii});
    B{ii} = MakeConvOperator(nx,ny,nz,S1);
end

C=MakeIdentityOperator(nx,ny,nz);   % it works better with the identity operator ???


% Compound 2
for ii=1:nacq
    S2 = MakeSpectrum(SpectrumC2{ii},directions{ii});
    A{ii} = MakeConvOperator(nx,ny,nz,S2);
end

% Construct the measurement operator
M = [];
for ii = 1:nacq
    M = [M;F{ii}*A{ii},F{ii}*B{ii}];
end


end


% Function definitions

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
