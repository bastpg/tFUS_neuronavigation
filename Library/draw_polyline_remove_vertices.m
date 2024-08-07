

function S = draw_polyline_remove_vertices( S )


debug = 0;

S2 = S;   S2.vertices(:,2) = [];

FIG = figure; 
patch( 'faces' , S2.faces , 'vertices' , S2.vertices , 'FaceColor' , [0.9 1 0.99] , 'FaceAlpha' , 1 );
title( 'Draw line' ); camlight; axis image off; 
pause(0.1);

PLINE = drawpolyline;

vx = S.vertices(:,1);
vy = S.vertices(:,3);

nsegs = length( PLINE.Position ) - 1;  % number of segments of the polyline

if debug == 1
    FIG1 = figure;
end

for ii = 1 : nsegs

    P1 = PLINE.Position(ii   , :);
    P2 = PLINE.Position(ii+1 , :);

    % P1 = [P1(2)  P1(1)];
    % P2 = [P2(2)  P2(1)];

    vylin = (vx-P1(1)) / (P2(1)-P1(1)) * (P2(2)-P1(2))  +   P1(2); 

    ind = find( vx>=(min(P1(1),P2(1)))  &  vx<=(max(P1(1),P2(1)))  &  vy<vylin  );
    S.faces = remove_faces_surf_mesh( S.faces , ind );

    if debug
        figure(FIG1);  clf(FIG1);
        patch( 'faces' , S.faces , 'vertices' , S.vertices , 'FaceColor' , [0.9 1 0.99] , 'FaceAlpha' , 1 );
        view([20 -10]); camlight; axis image off; hold on; pause(1);
    end

end


% only keep the largest surface
tmp = finddisconnsurf( S.faces );
nmax = -1;  iimax = -1;
for ii = 1 : numel(tmp)
    if size(tmp{ii},1) > nmax
        nmax = size(tmp{ii},1);
        iimax = ii;
    end
end
S.faces = tmp{iimax};

if debug
    figure(FIG1);  clf(FIG1);
    patch( 'faces' , S.faces , 'vertices' , S.vertices , 'FaceColor' , [0.9 1 0.99] , 'FaceAlpha' , 1 );
    view([20 -10]); camlight; axis image off; hold on; pause(1);
end















