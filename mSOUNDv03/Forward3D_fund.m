function [P_fundamental] = Forward3D_fund(mgrid, medium, excit_p, omega_c, ...
                           reflection_order, c_ref, varargin)
                       
% DESCRIPTION: Projecting the 3D wave field at the frequency of interest in
% the forward direction. Three different algorithms are considered here,
% which can be selected using the OPTIONAL INPUTS. For weakly inhomogeneous
% media (maximum sound velocity contrast less than 1.1, e.g., soft tissue),
% a mixed-domain method reported in IEEE UFFC 65(7), 1258–1267 (2018) is
% used. This is the default algorithm used in Forward3D_fund (solver 1).
% For strongly inhomogeneous media (e.g., skull), a split-step Fourier
% method with interpolation is used (solver 2). Note that this method can
% be also applied to weakly inhomogeneous media. Finally, for layered
% media (parallel layers of different acoustic media) regardless of the
% level of inhomogeneity, the conventional angular spectrum approach with
% analytically computed transmission and reflection coefficients is used
% (solver 3).

% USAGE:
% [P_fundamental] = Forward3D_fund(mgrid, medium, excit_p, omega_c)
% [P_fundamental] = Forward3D_fund(mgrid, medium, excit_p, omega_c, reflection_order)
% [P_fundamental] = Forward3D_fund(mgrid, medium, excit_p, omega_c, reflection_order, c_ref, ...)

% INPUTS:
% mgrid              Input structure to define the computational domain
% medium             Medium properties
% excit_p            Pressure distribution on the intial plane (excitation)
% omega_c            Center frequency 2*pi*fc
% reflection_order   The maximum order of reflection included in the simulation
% c_ref              A set of reference sound velocities covering the sound 
%                    velocity distribution of the acoustic medium.   

% Note: reflection_order is automatically set to 0 if there are only four
% inputs. Please only use this for homogeneous media, or weakly
% inhomogeneous media where you do not intend to compute the reflection.

% OPTIONAL INPUTS:
% 'correction'       Use the wave solver for strongly inhomogeneous media (solver 2) 
% 'layer'            Use the wave solver for layered media (solver 3)
% 'riemann'          Use a Riemann sum integration scheme for 
%                    a faster speed at the expense of a lower accuracy. 
%                    The default integation scheme is the trapezoidal rule. 
% 'NRL'              Activate the non-reflecting layer to reduce the
%                    wrap-around (spatial aliasing) error. 
% 'single'           Cast data type to single to reduce data size and possibly 
%                    speed up the computation at the expense of a loss in precision.
%
% Note: the specific order for the optional inputs is not important. For
% example, Forward3D_fund(mgrid, medium, excit_p, omega_c,
% reflection_order, c_ref,'correction','NRL','single') is the same as
% Forward3D_fund(mgrid, medium, excit_p, omega_c, reflection_order,
% c_ref,'single','correction','NRL')
%
% OUTPUTS:
% P_fundamental      Pressure distribution at the frequency of interest throughout the entire 3D domain

%% Check input 
if nargin<5
    reflection_order = 0;% no reflection is considered if input number is fewer than 5
    c_ref = [];
elseif nargin<6
    c_ref = [];
else 
    if ~isnumeric(c_ref)
error('The reference sound velocities c_ref cannot be found.')       
    end
end

MDMC = 0;  % flag for acoustic medium type
MDMA = 0;  % flag for non-reflecting layer
MDMI = 0;  % flag for integration scheme 
MDMD = 0;  % flag for data type (0: double; 1: single)
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;   % strongly inhomogeneous media
            case 'layer'
            MDMC = 2;   % layered media
            case 'NRL'
            MDMA = 1;   % non-reflecting layer
            case 'riemann' 
            MDMI = 1;
            case 'single'
            MDMD = 1;
            otherwise 
warning('At least one optional input is not recoganized. Please check the spelling.')                
        end
    end
end

%determine the level of medium inhomogeneity
if MDMC ==0&& max(max(max(medium.c)))/min(min(min(medium.c)))>1.05
warning('This may not be a weakly inhomogeneous medium. Please consider using the optional input "correction" or "layer", whichever is applicable.')
end

%determine if the medium is layered
if MDMC ==2
    if sum(sum(sum(abs(diff(medium.c,1,1)))))||sum(sum(sum(abs(diff(medium.c,1,2)))))~=0
warning('This might not be a layered medium. Please reconsider using optional input ''layer''. ')
    end 
end

if MDMA ==1
    M_gamma = repmat(1./(cosh(medium.NRL_alpha.*mgrid.abx_vec.')).^2, 1, mgrid.num_y) + ...
              repmat(1./(cosh(medium.NRL_alpha.*mgrid.aby_vec)).^2, mgrid.num_x, 1);
    M_gamma = M_gamma.*medium.NRL_gamma.*1i*omega_c;
    
%     figure; imagesc( abs(M_gamma) );  axis image off; colorbar; title('NRL');
%     pause(0.1);
    
else 
    M_gamma = 0;
end 

%Filter cutoff for removing singularities at points on the radiation circle
threshold = 10; % 10;% 10 is an empirical number. This number can be adjusted if needed. 

% preallocate an array to store reflected wave field
if reflection_order~=0
if MDMC~=2
p_ref = zeros(mgrid.num_x, mgrid.num_y,mgrid.num_z); 
elseif MDMC==2 
%identify layer interface 
layer_index1 = find(diff(medium.c(round(mgrid.num_x/2),round(mgrid.num_y/2),:)))+1;
layer_index2 = find(diff(medium.rho(round(mgrid.num_x/2),round(mgrid.num_y/2),:)))+1;
layer_index = unique(sort([layer_index1,layer_index2]));
p_ref = zeros(mgrid.num_x,mgrid.num_y,length(layer_index),'like',1+1i);
end
end

flag_homo = 0;     %flag for homogeneous media
if length(medium.rho) == 1&&length(medium.c) == 1
    flag_homo = 1; %identified as a homogeneous medium
end

if length(medium.rho) == 1 && length(medium.c) ~= 1
    medium.rho = medium.rho*ones(mgrid.num_x,mgrid.num_y,mgrid.num_z);
elseif length(medium.rho) ~= 1 && length(medium.c) == 1
    medium.c = medium.c*ones(mgrid.num_x,mgrid.num_y,mgrid.num_z);
end

%close all progress bars
% multiWaitbar( 'CloseAll' );
%%
if MDMC == 0
    % calculate the M term on the RHS of ODE
    [M_linear, Ktemp] = Mterm3D_fund(mgrid, medium, omega_c);
    
    % non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,1,mgrid.num_z);
     
    K = sqrt(Ktemp); %kz=sqrt(k^2-kx^2-ky^2)
    K(imag(K)==0)=-K(imag(K)==0);%K(imag(K)==0) is the propagating wave;K(imag(K)~=0) is the evanscent wave 
    
    % preallocate space to store the pressure at the fundamental frequency
    if MDMD == 1
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',single(1+1i));
    else
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',1+1i);    
    end
    
    %cast data type to "single"
    if MDMD == 1
    excit_p = single(excit_p);
    K = single(K);
    medium.c = single(medium.c);
    medium.rho = single(medium.rho);
    M_linear = single(M_linear);
    end
    
    % exponential terms to be used in the for loop
    expn2 = exp(1i*K*mgrid.dz); 
    expn3 = mgrid.dz*expn2./(2i.*K);
    
    f = excit_p;
    
    % Fourier transform of the input pressure W.R.T. x and y
    excit_F = fftshift(fft2(f)); 

    % multiWaitbar('Forward projection, please wait...',0);
    tic;
    for I = 2:mgrid.num_z
              
        M  = fftshift(fft2(M_linear(:,:,I-1).*f));       
        if MDMI ==1 %left-point Riemann sum scheme
        F1 = excit_F.*expn2 + expn3.*M;
        F1(isnan(F1)) = 0;
        F1(abs(K)<threshold )=0;%removing singularities at points on the radiation circle
        f = ifft2(ifftshift(F1));
        
        elseif MDMI ==0 %trapzoidal rule
        F1 = excit_F.*expn2 + expn3.*M;
        F1(isnan(F1)) = 0;
        F1(abs(K)<threshold )=0;
        
        f1  = ifft2(ifftshift(F1));      
        M1 = fftshift(fft2(M_linear(:,:,I).*f1));
        F2 = excit_F.*expn2 + 0.5*expn3.*(M + M1./expn2);
        F2(isnan(F2)) = 0;
        F2(abs(K)<threshold )=0;
        f = ifft2(ifftshift(F2));        
        end
                          
%% Computing the transmission coefficient (plane wave and layered-medium approximation)
if ~flag_homo
        T_rhoc = 2* medium.rho(:,:,I).*medium.c(:,:,I)./(medium.rho(:,:,I-1).*medium.c(:,:,I-1) + ...
                medium.rho(:,:,I).*medium.c(:,:,I));  
elseif flag_homo
    T_rhoc = 1;
end
%% Add reflection
        if reflection_order~=0&&any(T_rhoc-1,'all')
            p_ref(:,:, I-1) = f.*(T_rhoc-1); 
        end
        f = f.*T_rhoc;
        
        excit_F = fftshift(fft2(f));
        P_fundamental(:,:, I) = f;
        % multiWaitbar('Forward projection, please wait...','value',I/mgrid.num_z);
    end

elseif MDMC == 1

    c_ref = sort(c_ref);
    if length(c_ref)<2
    error('The number of reference sound velocities should be larger than one. For homogeneous media, the optional input "correction" should be disabled.');
    end
                       
    % calculate the M term on the RHS of ODE
    [M_linear, Ktemp] = Mterm3D_Mfund(mgrid, medium, omega_c, c_ref); 
    
    M_linear = M_linear + repmat(M_gamma,1,1,mgrid.num_z);
 
    K = sqrt(Ktemp); %kz=sqrt(k^2-kx^2-ky^2)
    K(imag(K)==0)=-K(imag(K)==0); %K(imag(K)==0) is the propagating wave;K(imag(K)~=0) is the evanscent wave
    
    % preallocate an array to store the pressure at the fundamental frequency
    if MDMD == 1
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',single(1+1i));
    else
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',1+1i);    
    end
    
    %interpolation factors
    H = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z,length(c_ref));
    
    if length(c_ref)==2
            Htemp = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp(medium.c>=c_ref(1)&medium.c<c_ref(2)) = (c_ref(2)-medium.c(medium.c>=c_ref(1)&medium.c<c_ref(2)))/(c_ref(2)-c_ref(1));
            H(:,:,:,1)=Htemp; 
            Htemp = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp(medium.c>=c_ref(1)&medium.c<=c_ref(2)) = (medium.c(medium.c>=c_ref(1)&medium.c<=c_ref(2))-c_ref(1))/(c_ref(2)-c_ref(1));
            H(:,:,:,2)=Htemp; 
    else 
         for II = 1:length(c_ref)
            if II==1
            Htemp = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp(medium.c>=c_ref(1)&medium.c<c_ref(2)) = (c_ref(2)-medium.c(medium.c>=c_ref(1)&medium.c<c_ref(2)))/(c_ref(2)-c_ref(1));
            H(:,:,:,II)=Htemp;
            elseif II>1&&II<length(c_ref)
            Htemp = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp(medium.c>=c_ref(II-1)&medium.c<c_ref(II)) = (medium.c(medium.c>=c_ref(II-1)&medium.c<c_ref(II))-c_ref(II-1))/(c_ref(II)-c_ref(II-1));
            
            Htemp1 = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp1(medium.c>=c_ref(II)&medium.c<c_ref(II+1)) = (c_ref(II+1)-medium.c(medium.c>=c_ref(II)&medium.c<c_ref(II+1)))/(c_ref(II+1)-c_ref(II));
            H(:,:,:,II)=Htemp+Htemp1; 
            else 
            Htemp = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z);
            Htemp(medium.c>=c_ref(II-1)&medium.c<=c_ref(II)) = (medium.c(medium.c>=c_ref(II-1)&medium.c<=c_ref(II))-c_ref(II-1))/(c_ref(II)-c_ref(II-1));
            H(:,:,:,II)=Htemp; 
            end
         end
    end
    
    %cast data type to "single"
    if MDMD == 1
    excit_p = single(excit_p);
    K = single(K);
    medium.c = single(medium.c);
    medium.rho = single(medium.rho);
    M_linear = single(M_linear);
    H = single(H);
    end
    
    % Fourier transform of the input pressure W.R.T. x and y 
    excit_F = fftshift(fft2(excit_p)); 
    f = excit_p;
    
    % exponential terms to be used in the for loop
    expn2 = exp(1i*K*mgrid.dz); 
    expn3 = mgrid.dz*expn2./(2i.*K);
     
    % multiWaitbar('Forward projection, please wait...',0);
    tic;
    for I = 2:mgrid.num_z
        
%         M  = fftshift(fft2(M_linear(:,:,I-1).*f));
%         f = 0;
%         for II = 1:length(c_ref) %split-step using different reference velocities
%         if MDMI ==1 %left-point Riemann sum scheme
%         F1 = excit_F.*expn2(:,:,II) + expn3(:,:,II).*M;
%         F1(isnan(F1)) = 0;
%         F1(abs(K(:,:,II))<threshold )=0;
%         f1 = ifft2(ifftshift(F1));
%         
%         elseif MDMI ==0 %trapzoidal rule
%         F1 = excit_F.*expn2(:,:,II) + expn3(:,:,II).*M;
%         F1(isnan(F1)) = 0;
%         F1(abs(K(:,:,II))<threshold )=0;
%         
%         f1  = ifft2(ifftshift(F1));      
%         M1 = fftshift(fft2(M_linear(:,:,I).*f1));
%         F2 = excit_F.*expn2(:,:,II) + 0.5*expn3(:,:,II).*(M + M1./expn2(:,:,II));
%         F2(isnan(F2)) = 0;
%         F2(abs(K(:,:,II))<threshold )=0;
%         f1 = ifft2(ifftshift(F2));
%         end
%         
%         % apply phase shift
%         f1 = f1*exp(1i*omega_c/c_ref(II)*mgrid.dz).*exp(-1i*omega_c./medium.c(:,:,I-1)*mgrid.dz);
%         % interpolation and summation
%         f1 = f1.*H(:,:,I-1,II);      
%         f = f+f1;
%         end
% %% Computing the transmission coefficient (plane wave and layered-medium approximation)   
%         T_rhoc = 2* medium.rho(:,:,I).*medium.c(:,:,I)./(medium.rho(:,:,I-1).*medium.c(:,:,I-1) + ...
%                 medium.rho(:,:,I).*medium.c(:,:,I));
% %% Add reflection 
%         if reflection_order~=0&&any(T_rhoc-1,'all')
%             p_ref(:,:, I-1) = f.*(T_rhoc-1); 
%         end       
%         f = f.*T_rhoc;    
%         
%         excit_F = fftshift(fft2(f));        
%         P_fundamental(:,:, I) = f;


        % SMALL MEMORY CODE:
        M  = fftshift(fft2(M_linear(:,:,I-1).*f));
        f = 0;
        for II = 1:length(c_ref) %split-step using different reference velocities
            M1 = excit_F.*expn2(:,:,II) + expn3(:,:,II).*M;
            M1(isnan(M1)) = 0;
            M1(abs(K(:,:,II))<threshold )=0;
            
            M1 = ifft2(ifftshift(M1));  % read_datFile_andCompare( M1 , 'M1_A_cref1.dat' , 1 , 256 , 256 );            
            M1 = fftshift(fft2(M_linear(:,:,I).*M1)); % read_datFile_andCompare( M1 , 'M1_B_cref1.dat' , 1 , 256 , 256 );
            M1 = excit_F.*expn2(:,:,II) + 0.5*expn3(:,:,II).*(M + M1./expn2(:,:,II));
            M1(isnan(M1)) = 0;
            M1(abs(K(:,:,II))<threshold )=0;
            M1 = ifft2(ifftshift(M1));  % read_datFile_andCompare( M1 , 'M1_C_cref1.dat' , 1 , 256 , 256 );
        
            % apply phase shift
            M1 = M1 * exp(1i*omega_c/c_ref(II)*mgrid.dz) .* exp(-1i*omega_c./medium.c(:,:,I-1)*mgrid.dz); 
            % interpolation and summation
            M1 = M1 .* H(:,:,I-1,II);      
            f = f + M1;  % read_datFile_andCompare( f , 'f_cref1.dat' , 1 , 256 , 256 );
        end
%% Computing the transmission coefficient (plane wave and layered-medium approximation)   
        T_rhoc = 2* medium.rho(:,:,I).*medium.c(:,:,I)./(medium.rho(:,:,I-1).*medium.c(:,:,I-1) + ...  
                medium.rho(:,:,I).*medium.c(:,:,I));
%% Add reflection 
        if reflection_order~=0&&any(T_rhoc-1,'all')
            p_ref(:,:, I-1) = f.*(T_rhoc-1); 
        end       
        f = f.*T_rhoc;    
        
        excit_F = fftshift(fft2(f));    % read_datFile_andCompare( excit_F , 'excit_F.dat' , 1 , 256 , 256 );
        P_fundamental(:,:, I) = f;
        
        
        %%%%%% end small memory code 

        % multiWaitbar('Forward projection, please wait...','value',I/mgrid.num_z);
        
    end
    
    elseif MDMC == 2
    % calculate the M term on the RHS of ODE
    [M_linear, Ktemp2] = Mterm3D_Mfund_layer(mgrid, medium, omega_c);
    
    % add non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,1,mgrid.num_z);
 
    K2 = sqrt(Ktemp2);%kz=sqrt(k^2-kx^2-ky^2)
    K2(imag(K2)==0)=-K2(imag(K2)==0); %K(imag(K)==0) is the propagating wave;K(imag(K)~=0) is the evanscent wave

    % preallocate an array to store the pressure at the fundamental frequency
    if MDMD == 1
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',single(1+1i));
    else
    P_fundamental = zeros(mgrid.num_x,mgrid.num_y,mgrid.num_z,'like',1+1i);    
    end
    
    %cast data type to "single"
    if MDMD == 1
    excit_p = single(excit_p);
    K2 = single(K2);
    medium.c = single(medium.c);
    medium.rho = single(medium.rho);
    M_linear = single(M_linear);
    end
    
    f = excit_p;
    % exponential terms to be used in the for loop
    expn2 = exp(1i*K2*mgrid.dz); 
    expn3 = mgrid.dz*expn2./(2i.*K2);
 
    % Fourier transform of the input pressure W.R.T. x and y 
    excit_F = fftshift(fft2(f)); 
    
    l_index = 1;
    % multiWaitbar('Forward projection, please wait...',0);
    tic;
    for I = 2:mgrid.num_z
                             
        M  = fftshift(fft2(M_linear(:,:,I-1).*f));       
        if MDMI ==1 %left-point Riemann sum scheme
        F1 = excit_F.*expn2(:,:,I-1) + expn3(:,:,I-1).*M;
        F1(isnan(F1)) = 0;
        F1(abs(K2(:,:,I-1))<threshold )=0;%removing singularities at points on the radiation circle
        excit_F = F1; 
        elseif MDMI ==0 %trapzoidal rule
        F1 = excit_F.*expn2(:,:,I-1) + expn3(:,:,I-1).*M;
        F1(isnan(F1)) = 0;
        F1(abs(K2(:,:,I-1))<threshold )=0;
        
        f1  = ifft2(ifftshift(F1));      
        M1 = fftshift(fft2(M_linear(:,:,I).*f1));
        F2 = excit_F.*expn2(:,:,I-1) + 0.5*expn3(:,:,I-1).*(M + M1./expn2(:,:,I-1));
        F2(isnan(F2)) = 0;
        F2(abs(K2(:,:,I-1))<threshold )=0;
        excit_F = F2;
        end
        
        %% Reflection and transmission computed using analytical solutions 
         if ismember(I,layer_index)
         theta_in = acos(-K2(:,:,I-1)./(omega_c./medium.c(:,:,I-1))); %incident angle 
         theta_tr = asin(medium.c(:,:,I)./medium.c(:,:,I-1).*sin(theta_in));%transmission angle-Snell's Law
         T = 2*medium.c(:,:,I).*medium.rho(:,:,I).*cos(theta_in)./(medium.c(:,:,I).*medium.rho(:,:,I).*...
         cos(theta_in)+medium.c(:,:,I-1).*medium.rho(:,:,I-1).*cos(theta_tr)); %transmission coefficient
         T(floor(real(theta_in)/pi*180)==90)=0;%remove singularities
         if reflection_order~=0
         p_ref(:,:, l_index) = ifft2(ifftshift(excit_F.*(T-1))); %reflected wave field
         end
         excit_F = excit_F.*T; 
         l_index = l_index+1;
         end 
         
        f = ifft2(ifftshift(excit_F)); 
        P_fundamental(:,:, I) = f;
        % multiWaitbar('Forward projection, please wait...','value',I/mgrid.num_z);
    end

end

%assigning intial pressure field to the 1st plane
P_fundamental(:,:,1) = excit_p;

%determine if reflection should be computed
if (max(max(max(medium.c)))    ~= min(min(min(medium.c)))) || ...
    (max(max(max(medium.rho))) ~= min(min(min(medium.rho))))
    reflection_option = 1;
else
    reflection_option = 0;
end

%run the reflection code
if reflection_order ~=0 && (reflection_option == 1)
    if MDMC == 0
        reflection = Reflection3D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, reflection_order, MDMI, MDMD);       
        P_fundamental = P_fundamental + reflection; 
    elseif MDMC == 1
        reflection = MReflection3D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, omega_c, reflection_order, c_ref, H, MDMI, MDMD); 
        P_fundamental = P_fundamental + reflection;        
    elseif MDMC == 2
        reflection = MReflection3D_fund_layer(mgrid, medium, p_ref, ...
                     M_linear, K2, expn2, expn3, omega_c, reflection_order, MDMI, MDMD);
        P_fundamental = P_fundamental + reflection;      
    end
end
disp(['mSOUND computation completed in ' num2str(toc) 's']); 
end