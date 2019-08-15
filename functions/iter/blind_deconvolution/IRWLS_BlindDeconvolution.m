function [k_out] = IRWLS_BlindDeconvolution(k_init, X, y, opts)

% solve % min 1/2\|Xk - Y\|^2 + \lambda \|k\|_1
% based on Dilip Krishnan Code 

if nargin == 3
  opts.lambda = 0;  
  % PCG parameters
  opts.pcg_tol = 1e-8;
  opts.pcg_its = 100;
  fprintf('Input options not defined - really no reg/constraints on the kernel?\n');
end

lambda = opts.lambda;
pcg_tol = opts.pcg_tol;
pcg_its = opts.pcg_its;
  

% precompute RHS
rhs=X'*y; 


% Set exponent for regularization
exp_a = 1;

k_out = k_init;
% outer loop
for iter = 1 : 1
  k_prev = k_out;
  % compute diagonal weights for IRLS
  weights_l1 = lambda .* (max(abs(k_prev), 0.0001) .^ (exp_a - 2));
  k_out = local_cg(k_prev, X,weights_l1, rhs, pcg_tol, pcg_its);
end
  
%==========================================================================

% local implementation of CG to solve the reweighted least squares problem
function k = local_cg(k, X, weights_l1, rhs, tol, max_its)

Ak = pcg_kernel_core_irls_conv(k, X,weights_l1);

r = rhs - Ak;

for iter = 1:max_its
  rho = (r(:)' * r(:));

  if (iter > 1)                      
    beta = rho / rho_1;
    p = r + beta*p;
  else
    p = r;
  end

  Ap = pcg_kernel_core_irls_conv(p, X, weights_l1);

  q = Ap;
  alpha = rho / (p(:)' * q(:) );
  k = k + alpha * p;                    
  r = r - alpha*q;                  
  rho_1 = rho;
  
  if (rho < tol)
    break;
  end;
end;

%==========================================================================

function out = pcg_kernel_core_irls_conv(k, X,weights_l1)
% This function applies the left hand side of the IRLS system to the
  
% first term: X'*X*k (quadratic term)
out_l2 = X'*X*k;
  
% second term: L1 regularization
out_l1 = weights_l1 .* k;
  
out = out_l2 + out_l1;



    