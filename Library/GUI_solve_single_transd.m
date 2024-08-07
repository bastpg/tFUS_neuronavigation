

function  [I_INT  PP  poro_int  c  rho  alpha ] = GUI_solve_single_transd( CENTER , NOR , x , y , z , poro , mgrid , medium , source_p , FREQ , reflection_order , varargin )
% function  [I_INT  PP  poro_int  c  rho  alpha ] = GUI_solve_single_transd( CENTER , NOR , x , y , z , poro , mgrid , medium , source_p , FREQ , reflection_order , varargin )
% 
% Computes beam intensity for a single transducer position using mSOUND
% IMPORTANT: It is imperative that the normal points toward the brain!!!

I_INT = [];
PP = [];
poro_int = [];
c   = [];
rho = [];
alpha = [];


c_ = CENTER;
n_ = NOR;

display = 0;
if numel(varargin) > 0
    S2 = varargin{1};
    display = 1;
end

if display == 1

    % plot skin volume, current point, transducer and transducer planes
    FIG = figure;
    P = patch( 'faces' , S2.faces , 'vertices' , S2.vertices , 'EdgeColor' , 'none' , 'FaceColor' , [0.2 0.2 0.2] , 'FaceAlpha' , 0.1 );
    set( FIG , 'Clipping' , 'off' ); camlight; axis image off; hold on;
    xlabel( 'X' ); ylabel('Y');  zlabel('Z');
    % quiver3(CENTER(:,1),CENTER(:,2),CENTER(:,3),NOR(:,1),NOR(:,2),NOR(:,3),'r');
    quiver3( c_(1) , c_(2) , c_(3) , n_(1) , n_(2) , n_(3) , 0.05 , 'g' , 'LineWidth' , 4 );
    plot3( c_(1) , c_(2) , c_(3) , 'og' , 'MarkerSize' , 10 , 'MarkerFaceColor' , 'green' );  
    U = cross( n_ , [ 1 1 1 ] );  U = U / norm(U);
    V = cross( n_ , U );          V = V / norm(V);
    theta = linspace( 0 , 2*pi , 100 )';
    xyzTrans = c_  +   60e-3/2 * cos(theta) * U  +  60e-3/2 * sin(theta) * V;
    plot3( xyzTrans(:,1) , xyzTrans(:,2) , xyzTrans(:,3) , '-g' , 'LineWidth' , 2 );    
    camzoom(1); pause(0.1); view( [ 180 0 ] );
end

% compute affine rotation matrix and interpolate scattering volume
R = compute_rot_matrix( n_ , [0 0 1] );  % the rotation matric aligns the transudcer normal onto [ 0;0;1 ]
T = [ 0 ; 0 ; 0] - R * c_';  %  the translation vector moveds the transducer center at [0;0;0];   

Raff = [ pinv(R)   -pinv(R)*T  ;
         0   0   0   1];

xf = mgrid.x;
yf = xf;
zf = mgrid.z(1:end-1);

poro_int = reshape( interpAffineGrid_MEX( x , y , z , permute( poro , [2 1 3] ) , xf , yf , zf , Raff ) , [ numel(xf)  numel(yf)  numel(zf) ] );

poro_int( find( poro_int<0 ) ) = 0;
poro_int( find( poro_int>1 ) ) = 1;

% Acoustics properties from CT units (IT'IS database)
c_brain = 1546;  
c_bone = 3514;  % m/s
rho_brain = 1046;  % Kg/m^3 (density)
rho_bone = 1908;  % kg/m^3 (density)    
alpha_brain= 0.5193;    % 6.8032  Np/m @ 1 MHz = 0.5193 dB/( MHz * cm )   using:   a = a/100 * 8.686 / (FREQ/1e6);
alpha_bone = 4.7385;    % 54.553  Np/m @ 1 MHz = 4.7385 dB/( MHz * cm )  

c =          c_brain * poro_int  +  c_bone     * ( 1 - poro_int );
rho =      rho_brain * poro_int  +  rho_bone   * ( 1 - poro_int );
alpha =  alpha_brain * poro_int  +  alpha_bone * ( 1 - poro_int );

if display == 1
    FIG2 = figure;
    subplot(2,2,1);  imagesc( squeeze(poro_int(:,round(end/2),:)) ); colorbar; axis image off; title('poro');
    subplot(2,2,2);  imagesc( squeeze(c(:,round(end/2),:)) ); colorbar; axis image off; title('c [m/s]');
    subplot(2,2,3);  imagesc( squeeze(alpha(:,round(end/2),:)) ); colorbar; axis image off; title('alpha [dB/(MHz*cm)]');
    subplot(2,2,4);  imagesc( squeeze(rho(:,round(end/2),:)) ); colorbar; axis image off; title('rho [kg/m^3]');
end


% run mSOUND simulation
medium.c( : , : , 2:end ) = c;
medium.rho( : , : , 2:end) = rho;
medium.ca( : , : , 2:end ) = alpha;

c_ref = [ min(medium.c(:))  max(medium.c(:)) ];

PP = Forward3D_fund( mgrid , medium , source_p , 2*pi*FREQ , reflection_order , c_ref , 'NRL' , 'correction' );

I = abs(PP).^2 ./ ( 2 * medium.rho .* medium.c ) * 1e6;  % Plane wave approximation

ind = find( medium.rho<1  |  medium.c<1 );
I(ind) = 0;

if display == 1
    figure(FIG2);
    subplot(2,2,4);  imagesc( squeeze(abs(I(:,round(end/2),:))) ); colorbar; axis image off; title('Pressure');
end

Raff_inv = [ R   T  ;
            0   0   0   1];

I_INT = reshape( interpAffineGrid_MEX( mgrid.x , mgrid.y , mgrid.z , I , x , y , z , Raff_inv ) , [ numel(x) numel(y)  numel(z) ] );
I_INT = permute( I_INT , [ 2 1 3 ] );

Raff_all = Raff;

if display == 1

    % plot beam in the global frame -- METHOD 1
    tmp = isosurface( x , y , z , I_INT , 0.25 );  
    [ V2   F2 ] = remesher( tmp.vertices,  tmp.faces, 0.001 , 1);

    figure(FIG); hold on;
    P = patch( 'faces' , F2 , 'vertices' , V2 , 'EdgeColor' , 'none' , 'FaceColor' , 'yellow' , 'FaceAlpha' , 1.0 );

    % % METHOD 2
    % tmp = isosurface( mgrid.x , mgrid.y , mgrid.z , abs(PP) , 0.50 );  
    % [ V2   F2 ] = remesher( tmp.vertices,  tmp.faces, 0.001 , 1);
    % tmp = [ V2 ones(size(V2,1),1) ] * Raff';  V2_ = tmp(:,1:3);
    % 
    % figure(FIG);
    % P = patch( 'faces' , F2 , 'vertices' , V2_ , 'EdgeColor' , 'none' , 'FaceColor' , 'blue' , 'FaceAlpha' , 0.5 );


end










