% Art Petrenko
% apetrenko@eos.ubc.ca
% March 2012
%
% For details see Donoho, Maleki and Montanari, "Message-passing algorithms
% for compressed sensing", 2009.

%% Definitions

maxiter = 1000;
tol = 1e-6;
options_spgl1.verbosity = 0;
options_spgl1.iterations = maxiter;

% Define an N-dimensional k-sparse signal, a measurement matrix and the
% measured data
N = 1000;
n = 100;
k = 10;
% undersampling
delta = n/N;
% sparsity
rho = k/n;

% Generate random sparse vector
x0 = zeros(N,1);
indices = randperm(N);
x0(indices(1:k)) = randn(k,1);

% mode 2 creates normalized columns in the Gaussian operator
A = opGaussian(n,N,2);
b = A*x0;

%% Recovery using different methods

% threshold parameter
lambda = 1;

% Solve
x_adj = A'*b;
[x_spgl1,~,~,info_spgl1] = spgl1(A,b,0,tol,[],options_spgl1);
% lsqr uses the relative residual norm as a stopping criterion
[x_lsqr,~,~,info_lsqr.iter,info_lsqr.r] = lsqr(A,b,tol/norm(b),maxiter);
[x_ist,info_ist] = ist(A,b,lambda,tol,maxiter);
[x_amp,info_amp] = ist(A,b,lambda,tol,maxiter,'amp');

% MSE
mse_adj = mse(x0,x_adj);
mse_spgl1 = mse(x0,x_spgl1);
mse_lsqr = mse(x0,x_lsqr);
mse_ist = mse(x0,x_ist);
mse_amp = mse(x0,x_amp);

% Succesful recovery criterion
suc_adj = DMMsuccess(x0,x_adj);
suc_spgl1 = DMMsuccess(x0,x_spgl1);
suc_lsqr = DMMsuccess(x0,x_lsqr);
suc_ist = DMMsuccess(x0,x_ist);
suc_amp = DMMsuccess(x0,x_amp);

% Plotting results
figure(1)
clf;
set(1,'Name','Sparse Recovery');
subplot(2,2,1)
plot(1:N,x_lsqr,'k', indices(1:k),x0(indices(1:k)),'ro');
title('LSQR Recovery')

subplot(2,2,2)
plot(1:N,x_spgl1,'k', indices(1:k),x0(indices(1:k)),'ro');
title('SPGl1 Recovery')

subplot(2,2,3)
plot(1:N,x_ist,'k', indices(1:k),x0(indices(1:k)),'ro');
title('IST Recovery')

subplot(2,2,4)
plot(1:N,x_amp,'k', indices(1:k),x0(indices(1:k)),'ro');
title('AMP Recovery')

figure(2)
clf;
set(2,'Name','Norm of residual and iterate: IST vs. AMP');
semilogy(0:info_ist.iter,info_ist.r,'g-',0:info_ist.iter,info_ist.s,'g--',...
         0:info_amp.iter,info_amp.r,'r-',0:info_amp.iter,info_amp.s,'r--');
xlabel('Iteration');
title('Norm of residual and threshold');
legend('IST ||r||','IST threshold / \lambda', ...
       'AMP ||r||','AMP threshold / \lambda');
axis tight