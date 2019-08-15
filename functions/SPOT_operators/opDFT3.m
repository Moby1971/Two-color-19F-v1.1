classdef opDFT3 < opOrthogonal
%OPDFT2  Two-dimensional fast Fourier transform (DFT).
%
%   opDFT2(M,N) creates a two-dimensional normalized Fourier transform
%   operator for matrices of size M by N. Input and output of the
%   matrices is done in vectorized form.
%
%   opDFT2(M,N,CENTERED) just like opDFT2(M,N), but with components
%   shifted to have to zero-frequency component in the center of the
%   spectrum, if the CENTERED flag is set to true.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties ( Access = private )
      funHandle         % Multiplication function
   end % properties
   
   properties ( SetAccess = private, GetAccess = public )
      inputdims         % Dimensions of the input
      centered          % Flag if operator created with center flag
   end % properties
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - Public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opDFT2. Constructor.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opDFT3(m,n,o,centered)
         if nargin <4 || nargin > 4
            error('Invalid number of arguments to opDFT2.');
         end
         
         op = op@opOrthogonal('DFT3',m*n*o,m*n*o);
         op.centered  = centered;
         op.cflag     = true;
         op.inputdims = [m,n,o];
         
         % Initialize function handle
         if centered
            op.funHandle = @opDFT3d_centered_intrnl;
         else
            op.funHandle = @opDFT3d_intrnl;
         end
      end % function opDFT2
      
   end % methods - public
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
      function y = multiply(op,x,mode)
         y = op.funHandle(op,x,mode);
      end % function multiply
      
   end % methods - protected

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - private
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = private )
       
       function y = opDFT3d_intrnl(op,x,mode)
           m = op.inputdims(1);
           n = op.inputdims(2);
           o = op.inputdims(3);
           
           if mode == 1
               y = reshape( fftn(reshape(full(x),m,n,o)) / sqrt(m*n*o), m*n*o, 1);
           else
               y = reshape(ifftn(reshape(full(x),m,n,o)) * sqrt(m*n*o), m*n*o, 1);
           end
       end % function opDFT3d_intrnl
      
      % Three-dimensional DFT - Centered
      function y = opDFT3d_centered_intrnl(op,x,mode)
          m = op.inputdims(1);
          n = op.inputdims(2);
          o = op.inputdims(3);
          
          if mode == 1
              y = fftshift(fftn(reshape(full(x),m,n,o))) / sqrt(m*n*o);
              y = reshape(y,m*n*o,1);
          else
              y = ifftn(ifftshift(reshape(full(x),m,n,o))) * sqrt(m*n*o);
              y = reshape(y,m*n*o,1);
          end
      end 
      
   end % methods - private
   
end % classdef
