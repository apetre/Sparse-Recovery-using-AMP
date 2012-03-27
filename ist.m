function [x,info] = ist(A,b,lambda,tol,maxiter,options)
% ist preforms iterative soft thresholding to solve A*x = b for x.
% 
% INPUTS
% lambda    The factor multiplying the "threshold" s, however the latter is
%           defined.
% tol       The stopping criterion that the 2-norm of the residual must satisfy.
% maxiter   The maximum number of iterations.
%
% OUTPUTS
% info.s    The values of the threshold at each iteration.
% info.r    The norm of the residual at each iteration.
% info.iter The number of iterations preformed. Note that each iteration
%           uses one multiply of A with x, and one of A' with x.

% The tolerance below which an entry is judged to be zero for the l0 norm.
L0tol = 1e-10;
L0x = 0;
j = 0;
[n,N] = size(A);
x = zeros(N,1);
r = A*x - b;
norm_r = norm(r);
info.s = zeros(maxiter+1,1);
info.r = zeros(maxiter+1,1);
info.s(1) = 0;
info.r(1) = norm_r;
if ~exist('options','var')
    options = 'ist';
end

while(j < maxiter && norm_r > tol)
    if strcmp(options,'amp')
        L0x = norm(double(abs(x)>L0tol),1);
    end
    % To make sure you understand the iteration indices: 
    % x_{k+1} = A'*r_k + x_k
    x = A'*r + x;
    x_sorted = sort(x,'descend');
    s = x_sorted(n);
    x = SoftThresh(x,s*lambda);
    % r_{k+1} = b - A*x_{k+1} + [ optional AMP term: r_{k}/n*norm(x_{k},0) ]
    r = b - A*x + r/n*L0x;
    norm_r = norm(r);
    j = j+1;
    info.r(j+1) = norm_r;
    info.s(j+1) = s;
end
info.iter = j;
info.s = info.s(1:j+1);
info.r = info.r(1:j+1);

end
