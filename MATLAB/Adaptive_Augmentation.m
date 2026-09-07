function [u_augmentation,W_hat_dot,psi] = Adaptive_Augmentation( ...
        estimated_state,estimated_tracking_error,W_hat,params,enabled)
%ADAPTIVE_AUGMENTATION Add the matched nonlinear term used in Task 5.
%   psi = [sin(2*phi_hat); sin(2*theta_hat);
%          sin(2*psi_hat); 1].
% This basis supplements Phi(x_hat) with the missing model uncertainty.

    if nargin < 5
        enabled = true;
    end

    estimated_state = estimated_state(:);
    estimated_tracking_error = estimated_tracking_error(:);
    if numel(estimated_state) ~= 6 ...
            || numel(estimated_tracking_error) ~= 6
        error('Adaptive_Augmentation:StateDimension', ...
            'State and tracking-error inputs must be 6-by-1.');
    end
    if ~isequal(size(W_hat),[4,3])
        error('Adaptive_Augmentation:EstimateDimension', ...
            'W_hat must be 4-by-3.');
    end

    psi = [sin(2*estimated_state(1));
           sin(2*estimated_state(2));
           sin(2*estimated_state(3));
           1];

    if enabled
        adaptation_signal = estimated_tracking_error.' ...
            *params.P*params.B;
        u_augmentation = -W_hat.'*psi;
        leakage_scale = params.task5.augmentation_sigma ...
            *norm(estimated_tracking_error,2);
        W_hat_dot = params.task5.Gamma_augmentation ...
            *psi*adaptation_signal ...
            -leakage_scale ...
                *params.task5.Gamma_augmentation*W_hat;
    else
        u_augmentation = zeros(3,1);
        W_hat_dot = zeros(4,3);
    end
end
