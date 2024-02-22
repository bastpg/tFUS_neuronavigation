classdef set_grid < handle
    
    % INPUT
    % dt         time step (only used for transient simulations)
    % t_length   total computational time (only used for transient simulations)
    % dx         spatial step size in the x-direction
    % dy         spatial step size in the y-direction
    % dz         spatial step size in the z-direction
    % x_length   computation domain size in the x direction
    % y_length   computation domain size in the y direction
    % z_length   computation domain size in the z direction
    
    % OUTPUT
    % mgrid   structure will be used in the simulation
    
    % STATIC PROPERTIES:
    % .x      coordinates in the x direction
    % .y      coordinates in the y direction
    % .z      coordinates in the z direction
    % .t      time array
    % .kx     wavevector in the x direction
    % .ky     wavevector in the y direction
    % .w      angular freuqency
    % .num_x   number of grid points in the x direction
    % .num_y   number of grid points in the y direction
    % .num_z   number of grid points in the z direction
    % .num_t   number of time steps
    
    % USAGE
    % mgrid = set_grid(dt, t_length, dx, x_length); 
    % mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length); 
    % mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length, dz, z_length); 
    
    % FUNCTION DESCRIPTION
    properties (GetAccess = public, SetAccess = private)
        % grid size in x-direction [grid points]
        num_x = 0;
        % grid size in y-direction [grid points]
        num_y = 0;
        % grid size in z-direction [grid points]
        num_z = 0;
        % grid size in t-direction [grid points]
        num_t = 0;
        % grid point spacing in x-direction [m]
      	dx = 0;
        % grid point spacing in y-direction [m]
        dy = 0;
        % grid point spacing in z-direction [m]
        dz = 0;
        % time step [s]
        dt = 0;
        % length in x-direction
        x_length = 0;
        % length in y-direction
        y_length = 0;
        % length in z-direction
        z_length = 0;    
        % length in the time domain
        t_length = 0;
        % x in the x-direction [m]  
        x;
        % y in the y-direction [m]  
        y;
        % z in the z-direction [m]  
        z;
        % time series
        t;
        % angular freuqency
        w; 
        % Nx x 1 vector of wavenumber components in the x-direction [rad/m]
        kx = 0;
        % Ny x 1 vector of wavenumber components in the y-direction [rad/m]
        ky = 0;             
        % number of dimensions
        dim = 0;  
        % vector for absorption layer in the x-direction [grid point]
        abx_vec = 0;
        % vector for absorption in the y-direction [grid point]
        aby_vec = 0;
    end
    
       properties(Dependent = true, GetAccess = public, SetAccess = private)
        % Nx x 1 vector of the grid coordinates in the x-direction [m] 
        x_vec;
        % Nx x 1 vector of the grid coordinates in the x-direction [m] 
        y_vec;
        % Nx x 1 vector of the grid coordinates in the x-direction [m] 
        z_vec;
        % Nx x 1 vector of wavenumber components in the x-direction [rad/m]
        kx_vec = 0;
        % Ny x 1 vector of wavenumber components in the y-direction [rad/m]
        ky_vec = 0;
        % angular frequency
        w_vec = 0;
        % vector for absorption layer in the x-direction [grid point]
        abx = 0;
        % vector for absorption in the y-direction [grid point]
        aby = 0;
 
       end 
    
       methods
           function mgrid = set_grid(varargin)  
               
               mgrid.dt = varargin{1};
               mgrid.t_length = varargin{2};
               % define time array
               mgrid.t = 0:mgrid.dt:mgrid.t_length; 
               % calculate the total number of points in time
               mgrid.num_t = length(mgrid.t);
               % call the function to define angular frequency omega
               mgrid.w = mgrid.makeOmega(mgrid.dt, mgrid.num_t);
               
               % assign the input values to the grid object
               switch nargin
                   case 4
                       % 1D uniform grid
                       mgrid.dx   = varargin{3};
                       mgrid.x_length = varargin{4};
                       
                       % calculate the grid points in x direction
                       mgrid.num_x = round(mgrid.x_length/mgrid.dx);
                       % set the number of dimensions
                       mgrid.dim = 1;   
                       
                   case 6
                       % 2D uniform grid
                       mgrid.dx      = varargin{3};
                       mgrid.x_length = varargin{4};
                       mgrid.dy      = varargin{5};
                       mgrid.y_length = varargin{6};  

                       % calculate the grid points in x- and y- directions
                       mgrid.num_x = round(mgrid.x_length/mgrid.dx);
                       mgrid.num_y = round(mgrid.y_length/mgrid.dy);
                       
                       % set the number of dimensions
                       mgrid.dim = 2; 
                       
                   case 8
                       
                       % 3D uniform grid
                       mgrid.dx      = varargin{3};
                       mgrid.x_length = varargin{4};
                       mgrid.dy      = varargin{5};
                       mgrid.y_length = varargin{6};        
                       mgrid.dz      = varargin{7};
                       mgrid.z_length = varargin{8};  
        
                       % calculate the grid points in x-, y- and
                       % z-directions
                       mgrid.num_x = round(mgrid.x_length/mgrid.dx);
                       mgrid.num_y = round(mgrid.y_length/mgrid.dy);        
                       mgrid.num_z = round(mgrid.z_length/mgrid.dz)+1;  
                       
                       % set the number of dimensions
                       mgrid.dim = 3; 
                   otherwise
                       error('Incorrect number of input arguments');
               end
               
               switch mgrid.dim
                   case 1
                       [mgrid.x, ~, ~] = mgrid.makeDim(mgrid.dx, mgrid.num_x);
                   case 2
                       [mgrid.x, mgrid.kx, mgrid.abx_vec] = mgrid.makeDim(mgrid.dx, mgrid.num_x);
                       [mgrid.y, ~, ~] = mgrid.makeDim(mgrid.dy, mgrid.num_y);
                       if mod(mgrid.num_x,2)==0 
                       warning('The number of grid points in the lateral direction is an even number. It is recommended that this number be an odd number');
                       end
                   case 3
                       [mgrid.x, mgrid.kx, mgrid.abx_vec] = mgrid.makeDim(mgrid.dx, mgrid.num_x);
                       [mgrid.y, mgrid.ky, mgrid.aby_vec] = mgrid.makeDim(mgrid.dy, mgrid.num_y);
                       %[mgrid.z, ~, ~] = mgrid.makeDim(mgrid.dz, mgrid.num_z);
                       mgrid.z = 0:mgrid.dz:(mgrid.num_z-1)*mgrid.dz;
                       if mod(mgrid.num_x,2)==0 || mod(mgrid.num_y,2)==0
                       warning('The number of grid points in the lateral direction is an even number. It is recommended that this number be an odd number');
                       end
               end  
           end 
       end

           
       % functions that can only be accessed by class members
       methods (Access = 'protected', Static = true) 
           
           % define the wavevector in the x/y/z direction: see doi: 10.1121/1.2812579
           function [x_vec, kx_vec, abx] = makeDim(dx, num_x)
               if mod(num_x,2) % odd num_x
                   x_vec   = (-num_x/2+1-1/2)*dx:dx:(num_x/2-1/2)*dx;
                   kx_vec  = [(-num_x/2+1-1/2):1:(num_x/2-1/2)]*2*pi/(num_x*dx);
                   abx     = [1:ceil(num_x/2)-1  ceil(num_x/2):-1:1];
               else % even number
                   x_vec   = (-num_x/2)*dx:dx:(num_x/2-1)*dx;
                   kx_vec  = [(-num_x/2):1:(num_x/2-1)]*2*pi/(num_x*dx);
                   abx     = [1:(num_x/2) (num_x/2):-1:1];
%                    x_vec   = (-num_x/2+1)*dx:dx:(num_x/2)*dx;
%                    kx_vec  = [(-num_x/2+1):1:(num_x/2)]*2*pi/(num_x*dx);
%                    abx     = [1:(num_x/2) (num_x/2):-1:1];             
               end
           end
           
           % define the angular frequency: https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
           function w_vec = makeOmega(dt, num_t)
               if mod(num_t,2) % odd length(t)
                   w_vec = 1/dt *2*pi* (-num_t/2+1-0.5:num_t/2-0.5)/num_t; %Angular frequency
               else %even length(t)
                   w_vec = 1/dt *2*pi* (-num_t/2:num_t/2-1)/num_t;%Angular frequency
          %w_vec(1)=-w_vec(1);  %the first frequency is the Nyquist frequency; 
               end
               
           end
           
       end

end