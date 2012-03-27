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

% Define an N-dimensional k-sparse signal, a measurement matrix and the data.
N = 1000;
% undersampling
delta = 0.4;
% sparsity
rho = 0.1;
n = floor(N*delta);
k = floor(n*rho);

% Generate random sparse vector
x0 = zeros(N,1);
indices = randperm(N);
x0(indices(1:k)) = randn(k,1);

% mode 2 creates normalized columns in the Gaussian operator
A = opGaussian(n,N,2);
b = A*x0;

% threshold parameter
lambda = 2.5;

%% Recovery using different methods

% Solve
[x_spgl1,~,~,info_spgl1] = spgl1(A,b,0,tol,[],options_spgl1);
[x_ist,info_ist] = ist(A,b,lambda,tol,maxiter);
[x_amp,info_amp] = ist(A,b,lambda,tol,maxiter,'amp');

% MSE
%mse_spgl1 = mse(x0,x_spgl1);
%mse_ist = mse(x0,x_ist);
%mse_amp = mse(x0,x_amp);

% Plotting results
figure(1)
clf;
subplot(2,2,1)
plot(1:N,x_spgl1,'k', indices(1:k),x0(indices(1:k)),'ro');
title('SPGl1')

subplot(2,2,3)
plot(1:N,x_ist,'k', indices(1:k),x0(indices(1:k)),'ro');
title('IST')

subplot(2,2,4)
plot(1:N,x_amp,'k', indices(1:k),x0(indices(1:k)),'ro');
title('AMP')

figure(2)
clf;
set(2,'Name','Norm of residual and iterate: IST vs. AMP');
semilogy(0:info_ist.iter,info_ist.r,'g-',0:info_ist.iter,info_ist.s,'g--',...
         0:info_amp.iter,info_amp.r,'r-',0:info_amp.iter,info_amp.s,'r--');
xlabel('Iteration');
legend('IST ||r||','IST \sigma', ...
       'AMP ||r||','AMP \sigma', 'Location','SouthEast');
ylim([1e-7,1]);
