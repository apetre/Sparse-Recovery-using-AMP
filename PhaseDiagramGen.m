maxiter = 1000;
tol = 1e-6;
options_spgl1.verbosity = 0;
options_spgl1.iterations = maxiter;
% Number of Monte Carlo runs per point of delta,rho phase space
MCnum = 4;
delta = linspace(0.02,0.9,16);
rho = linspace(0,0.5,16);

N = 1000;
lambda = 2.5;

phaseSpace_spgl1 = zeros(length(delta),length(rho));
phaseSpace_ist   = zeros(length(delta),length(rho));
phaseSpace_amp   = zeros(length(delta),length(rho));

for j_delta = 1:length(delta)
    for j_rho = 1:length(rho)
        for j_MC = 1:MCnum
            n = floor(N*delta(j_delta));
            k = floor(rho(j_rho)*n);
            % Generate random sparse vector
            x0 = zeros(N,1);
            indices = randperm(N);
            x0(indices(1:k)) = randn(k,1);
            % Generate matrix and data
            A = opGaussian(n,N,2);
            b = A*x0;
            
            % Solve
            [x_spgl1,~,~,info_spgl1] = spgl1(A,b,0,tol,[],options_spgl1);
            [x_ist,info_ist] = ist(A,b,lambda,tol,maxiter);
            [x_amp,info_amp] = ist(A,b,lambda,tol,maxiter,'amp');
            
            phaseSpace_spgl1(j_delta,j_rho) = ...
                phaseSpace_spgl1(j_delta,j_rho) + norm(x0-x_spgl1)/norm(x0);
            phaseSpace_ist(j_delta,j_rho) = ...
                phaseSpace_ist(j_delta,j_rho) + norm(x0-x_ist)/norm(x0);
            phaseSpace_amp(j_delta,j_rho) = ...
                phaseSpace_amp(j_delta,j_rho) + norm(x0-x_amp)/norm(x0);
        end
        phaseSpace_spgl1(j_delta,j_rho) = phaseSpace_spgl1(j_delta,j_rho)/MCnum;
        phaseSpace_ist(j_delta,j_rho)   = phaseSpace_ist(j_delta,j_rho)/MCnum;
        phaseSpace_amp(j_delta,j_rho)   = phaseSpace_amp(j_delta,j_rho)/MCnum;
    end
end
