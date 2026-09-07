function [Kx_dot,Kr_dot,Theta_dot] = eMod_MRAC( ...
        raw_Kx_dot,raw_Kr_dot,raw_Theta_dot, ...
        Kx_hat,Kr_hat,Theta_hat,tracking_error,params)
%EMOD_MRAC Scale leakage by ||e||_2, for example
% Kx_hat_dot=raw_Kx_dot-sigma_e*||e||_2*Gamma_x*Kx_hat.

    sigma_e = params.task4.sigma_e;
    if ~(isscalar(sigma_e) && isfinite(sigma_e) && sigma_e > 0)
        error('eMod_MRAC:InvalidLeakage', ...
            'params.task4.sigma_e must be a positive finite scalar.');
    end

    leakage_scale = sigma_e*norm(tracking_error,2);
    Kx_dot = raw_Kx_dot ...
        -leakage_scale*params.Gamma_x*Kx_hat;
    Kr_dot = raw_Kr_dot ...
        -leakage_scale*params.Gamma_r*Kr_hat;
    Theta_dot = raw_Theta_dot ...
        -leakage_scale*params.Gamma_theta*Theta_hat;
end
