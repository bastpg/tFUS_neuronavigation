function [M_linear, Ktemp] = Mterm3D_Mfund(mgrid, medium, omega_c, c_ref)

% DESCRIPTION:
% Computing the inhomogeneous term (M) at the frequency of interest
% (omega_c). Refer to Eq. 7 in J. Acoust. Soc. Am. 147 (6), 4055-4068,
% 2020. The sound velocity inhomogeneity is not considered here as it is
% instead considered in the wave propagation code (Forward3D_fund). This
% function is used by solver 2 which is applicable to both weakly and
% strongly inhomogeneous media.

% USAGE:
% [M_linear, Ktemp] = Mterm3D_Mfund(mgrid, medium, omega_c, c_ref)

% INPUTS:
% mgrid        Input structure to define the computational domain
% medium       Medium properties
% omega_c      Center frequency 2*pi*fc
% c_ref        A set of reference sound velocities covering the sound 
%              velocity distribution of the acoustic medium.   

% OUTPUTS:
% M_linear     M term (see code description)
% Ktemp        k_z^2 where k_z is the z-component of the wave-vector. These
%              values are computed for all reference sound velocities. 


%% 
% transform unit from dB/cm is Np/m
% alpha_Np  = abs(medium.ca.*(omega_c/2/pi/1e6).^medium.cb)*100*0.1151; 
alpha_Np  = abs(medium.ca.*(omega_c/2/pi/1e6))*100*0.1151;   

% diffusivity in the Westervelt eqn: delta = 2*c^3*alpha/omega^2
% compute absorption part in the M term: M_linear
M_linear  = 2i*alpha_Np*omega_c./medium.c;
%this part comes from (2*alpha_Np.*(medium.c.^3)/(omega_c^2))./medium.c.^4*1i*omega_c.^3

Ktemp = zeros(mgrid.num_x,mgrid.num_y,length(c_ref));
for I=1:length(c_ref)
Ktemp(:,:,I)     = (omega_c^2./c_ref(I)^2) -...
            (mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
            (ones(mgrid.num_x, 1)*mgrid.ky).^2;  
end
                
end
