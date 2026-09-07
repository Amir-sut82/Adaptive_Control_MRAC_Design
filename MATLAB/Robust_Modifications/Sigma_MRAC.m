function [Kx_dot,Kr_dot,Theta_dot] = Sigma_MRAC( ...
        raw_Kx_dot,raw_Kr_dot,raw_Theta_dot, ...
        Kx_hat,Kr_hat,Theta_hat,params)
%SIGMA_MRAC Add constant leakage, for example
% Kx_hat_dot=raw_Kx_dot-sigma*Gamma_x*Kx_hat.

    sigma = params.task4.sigma;
    if ~(isscalar(sigma) && isfinite(sigma) && sigma > 0)
        error('Sigma_MRAC:InvalidLeakage', ...
            'params.task4.sigma must be a positive finite scalar.');
    end

    Kx_dot = raw_Kx_dot-sigma*params.Gamma_x*Kx_hat;
    Kr_dot = raw_Kr_dot-sigma*params.Gamma_r*Kr_hat;
    Theta_dot = raw_Theta_dot ...
        -sigma*params.Gamma_theta*Theta_hat;
end
