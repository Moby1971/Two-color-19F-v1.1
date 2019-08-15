classdef opConvolve3D < opSpot
%OPCONVOLVE   One and two dimensional convolution operator.
%
%   opConvolve(M,N,KERNEL,OFFSET,MODE) creates an operator for one or
%   two-dimensional convolution, depending on the size of the KERNEL,
%   and the matrix or vector (MxN) the convolution is applied to. The
%   convolution is one dimensional only if KERNEL is a column vector
%   and N=1, or KERNEL is a row vector and M=1. The OFFSET parameter
%   determines the center of the KERNEL and has a default value of
%   [1,1]. When the OFFSET lies outside the size of the KERNEL, the
%   KERNEL is embedded in a zero matrix/vector with appropriate
%   center. For one-dimensional convolution, KERNEL may be a
%   scalar. Specifying an offset that is not equal to one where the
%   corresponding size of the kernel does equal one leads to the
%   construction of a two-dimensional convolution operator. There are
%   three types of MODE:
% 
%   MODE = 'regular'   - convolve input with kernel;
%          'truncated' - convolve input with kernel, but keep only
%                        those MxN entries in the result that
%                        overlap with the input;
%          'cyclic'    - do cyclic convolution of the input with a
%                        kernel that is wrapped around as many
%                        times as needed.
%
%   The output of the convolution operator, like all other
%   operators, is in vector form.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       funHandle     % Multiplication function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opConvolve3D(m,n,o,kernel,offset,mode)
  
           % warning('OFFSET NOT YET SUPPORTED')
           
           
          if nargin < 3
             error('opConvolve requires at least three parameters.');
          end

          if nargin < 4, offset = [];  end;
          if nargin < 5, mode = 'regular'; end;

          switch lower(mode)
             case {'cyclic'}
                cyclic = true; truncated = false;
    
             case {'truncated'}
              error(' truncated not supported')
    
             case {'default','regular',''}
              error(' default/regular not supported')
    
             otherwise
                error('Mode parameter must be one of ''regular'', ''cyclic'', or ''truncated''.');
          end

          % Check if one-dimensional convolution applies
          if (n == 1 && size(kernel,2) == 1) || (m == 1 && size(kernel,1) == 1)
             convolveOneDim = true;
             if length(offset) >= 2
                if ((n == 1 && size(kernel,2) == 1 && offset(2) ~= 1) || ...
                    (m == 1 && size(kernel,1) == 1 && offset(1) ~= 1))
                   convolveOneDim = false;
                else
                   offset = offset(1) * offset(2);
                end
             end
          else
             convolveOneDim = false;
             if length(offset) < 3
                error('The offset parameter needs to contain at least three entries.');
             end
          end


          if convolveOneDim
              error(' 1D not supported')
   
          else
             % ========= 3-dimensional case =========

             % Get basic information
             if isempty(offset), offset = [1,1,1];  end;
             if length(offset) < 3
                error('Offset parameter for 3D convolution needs to contain three entries.');
             end
             offset = offset(1:3);
             k      = [size(kernel,1),size(kernel,2),size(kernel,3)];
             cflag  = ~isreal(kernel);

             if cyclic
                % ========= Cyclic =========
      
                % Ensure offset(1) lies between 1 and m
                offset(1) = rem(offset(1)-1,m)+1;
                if offset(1) <= 0, offset(1) = offset(1) + m; end;

                % Ensure offset(2) lies between 1 and n
                offset(2) = rem(offset(2)-1,n)+1;
                if offset(2) <= 0, offset(2) = offset(2) + n; end;

                % Ensure offset(3) lies between 1 and o
                offset(3) = rem(offset(3)-1,o)+1;
                if offset(3) <= 0, offset(3) = offset(3) + o; end;
                
                % Wrap around kernel and zero pad if needed
                newKernel = zeros(m,n,o);
                for i=0:ceil(k(1)/m)-1
                   idx1 = 1:min(m,k(1)-i*m);
                   for j=0:ceil(k(2)/n)-1
                      idx2 = 1:min(n,k(2)-j*n);
                      for q=0:ceil(k(3)/o)-1
                        idx3 = 1:min(o,k(3)-q*o);
                        newKernel(idx1,idx2,idx3) = newKernel(idx1,idx2,idx3) + kernel(i*m+idx1,j*n+idx2,q*o+idx3);
                      end
                   end
                end
                kernel = newKernel;
       
                %%%%%%%TO DO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply cyclic shifts to correct for offset
%                 kernel = [kernel(offset(1):end,offset(2):end), kernel(offset(1):end,1:offset(2)-1); ...
%                           kernel(1:offset(1)-1,offset(2):end), kernel(1:offset(1)-1,1:offset(2)-1)];
%                 kernel=circshift(kernel,offset(1)-1,1);
%                 kernel=circshift(kernel,offset(2)-1,2);
%                 kernel=circshift(kernel,offset(3)-1,3);

                
                % Precompute kernel in frequency domain
                fKernel = fftn(full(kernel));

                % Create function handle and determine operator size
                fun   = @(x,mode) opConvolveCircular3D_intrnl(fKernel,m,n,o,cflag,x,mode);
                nRows = m*n*o;
                nCols = m*n*o;
             else
              error('not supported')
             end
          end

          % Construct operator
          op = op@opSpot('Convolve3D', nRows, nCols);
          op.cflag     = cflag;
          op.funHandle = fun;
       end % Constructor

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = op.funHandle(x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%======================================================================

function y = opConvolveCircular3D_intrnl(fKernel,m,n,o,cflag,x,mode)
if mode == 1
   y = ifftn(fKernel.*fftn(full(reshape(x,m,n,o))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
else
   y = ifftn(conj(fKernel).*fftn(full(reshape(x,m,n,o))));
   y = y(:);
   if (~cflag && isreal(x)), y = real(y); end;
end
end
