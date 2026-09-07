function [Am,Bm,Cm] = Reference_Model(params)
%REFERENCE_MODEL Build xmdot=Am*xm+Bm*r and ym=Cm*xm.

    Am = [zeros(3),eye(3);
          -params.Kp,-params.Kd];

    Bm = [zeros(3);
          params.Kp];

    Cm = [eye(3),zeros(3)];
end
