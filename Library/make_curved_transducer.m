

function   TD_mask = make_curved_transducer( x , y , z , C , N , R , apDiam )
% function   TD_mask = make_curved_transducer( grid , C , N , R , apDiam )
% 
% x , y , z: grid
% C: center of transducer = backmost point of the curved surface
% N: transducer normal 
% R: transducer focal radius 
% apDiam: aperture diameter



addpath( '/space/guerin/CODE/iso2mesh/' );
% addpath( 'Z:\CODE\iso2mesh' );

N = N / norm(N);  

% center of the sphere
O = C + N * R;


[ xx  yy  zz ]  =  ndgrid( x , y , z );

temp = sqrt( ( xx - O(1) ).^2  +  ( yy - O(2) ).^2   +  ( zz - O(3) ).^2 );

mask1 = temp <= R;  % transducer sphere
mask2 = thinbinvol( mask1 , 1 );   % erode by 1 voxel to create shell

TD_mask = mask1 - mask2;

% remove voxels not in the aperture
DIST = [  xx(:) - C(1)    yy(:) - C(2)     zz(:) - C(3)  ] * N';

dAP = R  -  sqrt( R^2 - apDiam^2/4 );
TD_mask( DIST > dAP ) = 0;





    
    




