function [Kx_dot,Kr_dot,Theta_dot,is_active] = Deadzone_MRAC( ...
        raw_Kx_dot,raw_Kr_dot,raw_Theta_dot,tracking_error,epsilon)
%DEADZONE_MRAC Freeze every estimate while ||e||_2<=epsilon.

    if ~(isscalar(epsilon) && isfinite(epsilon) && epsilon >= 0)
        error('Deadzone_MRAC:InvalidThreshold', ...
            'epsilon must be a finite nonnegative scalar.');
    end

    is_active = norm(tracking_error,2) > epsilon;
    if is_active
        Kx_dot = raw_Kx_dot;
        Kr_dot = raw_Kr_dot;
        Theta_dot = raw_Theta_dot;
    else
        Kx_dot = zeros(size(raw_Kx_dot));
        Kr_dot = zeros(size(raw_Kr_dot));
        Theta_dot = zeros(size(raw_Theta_dot));
    end
end
