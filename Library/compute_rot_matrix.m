


function R = compute_rot_matrix( dir1 , dir2 )
% function R = compute_rot_matrix( dir1 , dir2 )
% compute the 3x3 ritation matrix that aligns dir1 onto dir2



% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d


dir1 = dir1'/norm(dir1);  % norm = 1 
dir2 = dir2'/norm(dir2);  

if norm( dir1 - dir2 ) < 1e-6  % formula below does not handle well the case dir1 == dir2
    R = eye(3);
else
    GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; ...
              norm(cross(A,B)) dot(A,B)  0; ...
              0              0           1];

    FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

    UU = @(A,B) ( FFi(A,B)*GG(A,B)*inv(FFi(A,B)) ) ;

    R = UU( dir1 , dir2 );  % R is the rotation that transforms dir2 (the current direction of the transducer) to dir1 (direction in the water simulation)
end

% verify calculation
if norm( R*dir1 -  dir2 ) > 1e-5
    error( 'Incorrect rotation calculation!' );
end

