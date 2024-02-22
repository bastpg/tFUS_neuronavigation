



function values_sm = smooth_faces( F , V , values )
% function values_sm = smooth_faces( F , V , values )
% average face values over faces connected once (1 vertex in common)

values_sm = zeros(size(values));

for ii = 1 : size(F,1)

    ind = find( F == F(ii,1) |  F == F(ii,2) |  F == F(ii,3) );
    [ face_num  iy ] = ind2sub( size(F) , ind );

    face_num = unique(face_num);  % remove duplicates

    face_pos  = ( V(F(ii,1),:) + V(F(ii,2),:) + V(F(ii,3),:) ) / 3;  % center of current face
    face_pos2 = ( V(F(face_num,1),:) + V(F(face_num,2),:) + V(F(face_num,3),:) ) / 3;  % center of triangles
    dist = sqrt( sum( abs( face_pos - face_pos2 ).^2  , 2 ) );

    weights = 1./dist;
    ind2 = find(dist>0);  % this is to remove the currenrt face from the list, which would have INF weight

    values_sm(ii) = 0.5 * values(ii)  +   0.5 * sum( values(face_num(ind2)) .* weights(ind2) / sum(weights(ind2)) );

end




































