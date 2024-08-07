

function vert_vals = interp_face2vert( faces ,  nodes , face_vals , list_nodes )
% function vert_vals = interp_face2vert( faces ,  nodes , face_vals , list_nodes )


nnodes = size(nodes,1);
nfaces = size(faces,1);

vert_vals = zeros(nnodes,1);

for ii = 1 : numel(list_nodes)

    if mod(ii,100)==0
        fprintf('[%d/%d]\n' , ii , numel(list_nodes))
    end
    IND = find( faces(:) == list_nodes(ii) );

    weights = [];
    vals = [];
    for jj = 1 : numel(IND)
        face_num = mod( IND(jj) , nfaces );
        if face_num == 0
            face_num = nfaces;
        end
        face_cent = ( nodes(faces(face_num,1),:) + nodes(faces(face_num,2),:) + nodes(faces(face_num,3),:) ) / 3;  % center of triangle
        dist2node = sqrt( sum( abs( face_cent - nodes(list_nodes(ii),:) ).^2 ) );
        weights = [ weights ; 1/dist2node ];
        vals = [ vals ;  face_vals(face_num) ];
    end

    vert_vals(list_nodes(ii)) = sum( vals .* weights ) / sum(weights);


end












