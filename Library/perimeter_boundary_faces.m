

function P = perimeter_boundary_faces( BE )

P= [];

while ~isempty(BE)

    pc = BE(1,2);  % we start at the first point of the first edge 
    pe = BE(1,1);  % the last point has to be the second point of the first edge (closed perimeter)
    % add first edge to perimeter and remove it from further consideration
    P{end+1} = BE(1,:);
    BE(1,:) = [];  
    
    while pc ~= pe  % add edges until the current point is different from the last
    
        % find edge sharing current point
        a = find( BE(:,1) == pc );
        b = find( BE(:,2) == pc );
        if      ~isempty(a) &   isempty(b)
            i = a;
            e2a = [ BE(i,1)  BE(i,2) ];
        elseif   isempty(a) &  ~isempty(b)
            i = b;
            e2a = [ BE(i,2)  BE(i,1) ];
        else
            error('Could not find valid next edge in perimeter_boundary_faces().');
        end
    
        % add this edge to perimeter and remove it from further consideration
        P{end} = [ P{end} ; e2a ];
        BE(i,:) = [];
    
        % update current point
        pc = e2a(2);
    
    end

end

fprintf( 'Found %d disconnected perimeters!\n' , numel(P) );


















