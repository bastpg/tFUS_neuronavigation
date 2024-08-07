

function imgout = reslice_( imgin )
% reslice volume for better visualization inside the navigation app

tmp = permute( imgin , [ 1 3 2] );
imgout = tmp( : , : , end:-1:1 );






