function [M_linear, Ktemp2] = Mterm3D_Mfund_layer(mgrid, medium, omega_c)

% DESCRIPTION:
% Computing the inhomogeneous term (M) at the frequency of interest
% (omega_c). Refer to Eq. 7 in J. Acoust. Soc. Am. 147 (6), 4055-4068,
% 2020. Note that the nonlinear term associated with beta is not considered
% here. The sound velocity inhomogeneity is also not considered here as it
% is instead considered in the wave propagation code (Forward3D_fund). This
% function is used by solver 3 which is applicable to layered media
% (parallel layers of different acoustic media) regardless of the level of
% inhomogeneity.

% USAGE:
% [M_linear, Ktemp2] = Mterm3D_Mfund_layer(mgrid, medium, omega_c)

% INPUTS:
% mgrid        Input structure to define the computational domain
% medium       Medium properties
% omega_c      Center frequency 2*pi*fc

% OUTPUTS:
% M_linear     M term
% Ktemp2       k_z^2 where k_z is the z-component of the wave-vector. These
%              values take into account the different sound velocities in
%              different layers. 

%% 
% transform unit from dB/cm is Np/m
alpha_Np  = abs(medium.ca.*(omega_c/2/pi/1e6).^medium.cb)*100*0.1151;  

% diffusivity in the Westervelt eqn: delta = 2*cw^3*alpha/w^2
% compute absorption part in the M term: M_linear
 M_linear  = 2i*alpha_Np*omega_c./medium.c;
%this part comes from (2*alpha_Np.*(medium.c.^3)/(omega_c^2))./medium.c.^4*1i*omega_c.^3

kxy = squeeze((mgrid.kx.'*ones(1, mgrid.num_y)).^2 +...
             (ones(mgrid.num_x, 1)*mgrid.ky).^2);
        
Ktemp2    = (omega_c^2./medium.c.^2) - repmat(kxy, 1,1,mgrid.num_z);  
end