


function P = rayleighSommerfeldIntegral( transdPoints , transdCAmp , x , y , z , c , rho , freq )
% function P = rayleighSommerfeldIntegral( transdPoints , transdCAmp , x , y , z , c , rho , freq )
% transdPoints: N1*3   list of points with a US emmitter ([ x y z ] for each point)
% transdCAmp: list of complex amplitudes for each US emitter point
% x , y , z: list of computation points x, y and z coordinate (have the same length)
% c, rho, freq: m/s, kg/m^3, Hz
% P (output): pressure at the computation points

omega = 2 * pi * freq;

% [ a  b ] = ndgrid( x , transdPoints(:,1) );
% DX = a - b;
% 
% [ a  b ] = ndgrid( y , transdPoints(:,2) );
% DY = a - b;
% 
% [ a  b ] = ndgrid( z , transdPoints(:,3) );
% DZ = a - b;
% 
% R = sqrt( DX.^2  +   DY.^2  +   DZ.^2 );
% 
% P = ( 1./R * 1j * rho * c * omega / (2 * pi) .*   exp( -1j * omega * R / c ) )   *    transdCAmp;


P = zeros( size(x) );

for ii = 1 :  size( transdPoints , 1 )  
    
    % distance between the transducer (current point) and the observation
    % points
    R = sqrt( (x-transdPoints(ii,1)).^2  +   (y-transdPoints(ii,2)).^2  +  (z-transdPoints(ii,3)).^2 );
            
    P  =   P   +   transdCAmp(ii) * 1j * rho * c * omega / (2 * pi) * exp( -1j * omega * R / c ) ./ R;
    P( find( R == 0 ) ) = 0;
    
end








