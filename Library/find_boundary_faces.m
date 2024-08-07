



function [ BF  BE ] = find_boundary_faces( F )
% BF: boundary faces
% BE: boundary edges


A1 = F(:,1);
A2 = F(:,2);
A3 = F(:,3);

BF = [];
BE = [];

for ii = 1 : length(F)

    % edges of current face
    e1 = [ F(ii,1)  F(ii,2) ];
    e2 = [ F(ii,1)  F(ii,3) ];
    e3 = [ F(ii,2)  F(ii,3) ];

    % find faces that share e1
    N = ( A1==e1(1) |  A1==e1(2) )  +  ( A2==e1(1) |  A2==e1(2) )  +  ( A3==e1(1) |  A3==e1(2) );
    indF1 = find(N>=2);  % the faces that share this edge share the two points of that edge

    % find faces that share e2
    N = ( A1==e2(1) |  A1==e2(2) )  +  ( A2==e2(1) |  A2==e2(2) )  +  ( A3==e2(1) |  A3==e2(2) );
    indF2 = find(N>=2); 

    % find faces that share e3
    N = ( A1==e3(1) |  A1==e3(2) )  +  ( A2==e3(1) |  A2==e3(2) )  +  ( A3==e3(1) |  A3==e3(2) );
    indF3 = find(N>=2);

    % boundary faces are those with less than 3 shared faces (+1 for the current face)
    indAll = unique( [ indF1 ; indF2 ; indF3 ] );
    if numel(indAll) < 4
        BF = [ BF ; ii ];
    end

    % boundary edges
    if numel(indF1) == 1
        BE = [ BE ; e1 ];
    end
    if numel(indF2) == 1
        BE = [ BE ; e2 ];
    end
    if numel(indF3) == 1
        BE = [ BE ; e3 ];
    end
 

end






















