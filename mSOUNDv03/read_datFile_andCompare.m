



function read_datFile_andCompare( M , fname , flag_complex , nx , ny )
% function read_datFile_andCompare( M , fname , flag_complex , nx , ny )


if flag_complex == 1
    
    M_ = reshape( read_ddata( fname , -1 ) , nx , ny , 2 );
    
    figure; 
    subplot( 2 , 3 , 1 ); imagesc( real(M) );  axis image off; colorbar;  title( 'REAL -- ref' )
    subplot( 2 , 3 , 2 ); imagesc( M_(:,:,1) );  axis image off; colorbar;  title( 'REAL -- mex' )
    subplot( 2 , 3 , 3 ); imagesc( real(M) - M_(:,:,1) );  axis image off; colorbar;  title( 'REAL -- (ref-mex)' )
    
    subplot( 2 , 3 , 4 ); imagesc( imag(M) );  axis image off; colorbar;  title( 'IMAG -- ref' )
    subplot( 2 , 3 , 5 ); imagesc( M_(:,:,2) );  axis image off; colorbar;  title( 'IMAG -- mex' )
    subplot( 2 , 3 , 6 ); imagesc( imag(M) - M_(:,:,2) );  axis image off; colorbar;  title( 'IMAG -- (ref-mex)' )

else
    
    M_ = reshape( read_ddata( fname , -1 ) , nx , ny );
    
    figure; 
    subplot( 2 , 3 , 1 ); imagesc( real(M) );  axis image off; colorbar;  title( 'ref' )
    subplot( 2 , 3 , 2 ); imagesc( M_ );  axis image off; colorbar;  title( 'mex' )
    subplot( 2 , 3 , 3 ); imagesc( real(M) - M_ );  axis image off; colorbar;  title( 'ref-mex' )
    

    
end
















