% Art Petrenko
% apetrenko@eos.ubc.ca
% March 2012
%
% For details see Donoho, Maleki and Montanari, "Message-passing algorithms
% for compressed sensing", 2009.

%% Definitions

L0tol = 1e-10;
maxiter = 1000;
tol = 1e-6;
options_spgl1.verbosity = 0;
options_spgl1.iterations = maxiter;

% Define an N-dimensional k-sparse signal, a measurement matrix and the
% measured data
N = 1000;
n = 100;
k = 5;
% undersampling
delta = n/N;
% sparsity
rho = k/n;

x0 = zeros(N,1);
indices = randperm(N);
x0(indices(1:k)) = randn(k,1);

% mode 2 creates normalized columns in the Gaussian operator
A = opGaussian(n,N,2);
b = A*x0;

%% Recovery using different methods

% threshold parameter
lambda = 3;
s = zeros(maxiter,1);
norm_r = zeros(maxiter,1);
norm_x = zeros(maxiter,1);

x_adj = A'*b;
x_spgl1 = spgl1(A,b,0,tol,[],options_spgl1);
x_lsqr = lsqr(A,b,tol,maxiter);

x = zeros(N,1);
r = b - A*x;
norm_r(1) = norm(r);
j = 2;
while(j < maxiter)
    L0 = norm(double(x>L0tol),1);
    x = A'*r + x;
    x_sorted = sort(x,'descend');
    s(j-1) = x_sorted(n);
    x = softThresh(x,lambda*s(j-1));
    r = b - A*x + r/n*L0;
    norm_r(j) = norm(r);
    norm_x(j) = norm(x);
    j = j+1;
end
x_ist = x;

% MSE
mse_adj = mse(x0,x_adj);
mse_spgl1 = mse(x0,x_spgl1);
mse_lsqr = mse(x0,x_lsqr);
mse_ist = mse(x0,x_ist);

% Plotting results
figure(1)
clf;
set(1,'Name','Sparse Recovery');
subplot(2,2,1)
plot(1:N,x_adj,'ko', indices(1:k),x0(indices(1:k)),'ro');
title('A^{H} b')

subplot(2,2,2)
plot(1:N,x_spgl1,'k', indices(1:k),x0(indices(1:k)),'ro');
title('SPGl1 Recovery')

subplot(2,2,3)
plot(1:N,x_lsqr,'k', indices(1:k),x0(indices(1:k)),'ro');
title('LSQR Recovery')

subplot(2,2,4)
plot(1:N,x,'k', indices(1:k),x0(indices(1:k)),'ro');
title('IST Recovery')

figure(2)
clf;
set(2,'Name','Norm of residual and iterate');
loglog(1:j-1,norm_r(1:j-1), 1:j-1,norm_x(1:j-1), 1:j-1,s(1:j-1));
legend('r','x', 's');
axis tight

% Raise first figure
figure(1)