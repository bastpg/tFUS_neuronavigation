


function   FF = remove_faces_surf_mesh( FF , ind )
% function   FF = remove_faces_surf_mesh( FF , ind )
% remove faces connected to vertices with indices "ind" 

F2R = [];

for ii = 1 : numel(ind)
   
    F1 = find(  FF(:,1) == ind(ii) );
    F2 = find(  FF(:,2) == ind(ii) );
    F3 = find(  FF(:,3) == ind(ii) );
    
    F2R = [ F2R  ;   unique([ F1  ; F2 ; F3 ]) ];
    
end

FF( unique(F2R) , : ) = [];







