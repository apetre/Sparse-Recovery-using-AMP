% Art Petrenko
% apetrenko@eos.ubc.ca
% March 2012

%% Sparse Recovery

max_iter = 100;
tol = 1e-6;
options_spgl1.verbosity = 0;
options_spgl1.iterations = max_iter;

% Define an N-dimensional k-sparse signal, a measurement matrix and the
% measured data
N = 1000;
n = 100;
k = 5;

x0 = zeros(N,1);
indices = randperm(N);
x0(indices(1:k)) = randn(k,1);

A = opGaussian(n,N);
b = A*x0;

% Recovery using different methods
x_adj = A'*b;
x_spgl1 = spgl1(A,b,0,tol,[],options_spgl1);
x_lsqr = lsqr(A,b,tol,max_iter);

% MSE
mse_adj = mse(x0,x_adj);
mse_spgl1 = mse(x0,x_spgl1);
mse_lsqr = mse(x0,x_lsqr);

%%
figure(1)
clf;
set(1,'Name','Sparse Recovery');
subplot(2,2,1)
plot(1:N,x_adj,'k', indices(1:k),x0(indices(1:k)),'ro');
title('A^{H} b')

subplot(2,2,2)
plot(1:N,x_spgl1,'k', indices(1:k),x0(indices(1:k)),'ro');
title('SPGl1 Recovery')

subplot(2,2,3)
plot(1:N,x_lsqr,'k', indices(1:k),x0(indices(1:k)),'ro');
title('LSQR Recovery')