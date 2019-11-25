function M = MakeMeasurementOperator2C(SpectrumC1,SpectrumC2,nx,ny,nz,directions,varargin)
% Constructs the Matrix (F 0/0 F)(C 0/0 C)

nacq = length(directions); % number of acquisitions in system

if nargin>6
    F = varargin{1};
else
    for ii = 1:nacq
        F{ii} = opDFT3(nx,ny,nz,1);
    end
end


% Compound 1
for ii=1:nacq
    
    S1 = MakeSpectrum(SpectrumC1,directions{ii});
    B{ii} = MakeConvOperator(nx,ny,nz,S1);
    
end

% C = MakeIdentityOperator(nx,ny,nz);   % used previously for PFCE


% Compound 2
for ii=1:nacq
    
    S2 = MakeSpectrum(SpectrumC2,directions{ii});
    A{ii} = MakeConvOperator(nx,ny,nz,S2);
    
end

% Construct the measurement operator
M = [];
for ii = 1:nacq
    
    M = [M;F{ii}*A{ii},F{ii}*B{ii}];
    
end


end


% Function definitions

% Deconvolution operator
function A = MakeConvOperator(nx,ny,nz,S)

A = opConvolve3D(nx,ny,nz,S,[1 1 1],'cyclic');

end


% Identity operator
function B = MakeIdentityOperator(nx,ny,nz)

B = opDirac(nx*ny*nz);

end


% Flipping the spectrum according to the measurement directions
function S = MakeSpectrum(Spectrum,dir)

switch dir
    
    case 'HF'
        S = Spectrum;
        
    case 'FH'
        S = flip(circshift(Spectrum,-1));
        
    case 'LR'
        S = Spectrum';
        
    case 'RL'
        S = flip(circshift(Spectrum,-1))';
        
    otherwise
        error('direction error');
        
end

end
