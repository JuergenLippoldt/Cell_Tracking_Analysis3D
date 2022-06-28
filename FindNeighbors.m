function [VertexNeighbors, CellNeighbors] = FindNeighbors( C )
% C is a cell with tesselation cells, i.e.
% C{1} conains which vertices cell 1 is made up of.
% return value: Neighbors{n} tells which cells are neighbors of cell n

	% first find out which vertex belongs to which cell
	VertexNeighbors = cell(max(cellfun(@max, C)), 1);
	for c = 1:length(C)
		for v = 1:length(C{c}) 
			VertexNeighbors{C{c}(v)}(end+1) = c;		
		end
	end
	% In VertexNeighbors, every vertex should now have three cell indices
	% it belongs to
    VertexNeighbors{1}=[];  %% first Vertex is Inf - gives wrong neighbors

	% Now for each cell, find the unique neighbors of all vertices the cell
	% is made up of
	CellNeighbors = cellfun( @(x) unique([VertexNeighbors{x}]), C, 'uniformoutput', false );
	% be aware that the cell index of the cell itself is contained in that
	% info!
    
    
    %% I still want to know how boarders Inf - to sort them out
    VertexNeighbors = cell(max(cellfun(@max, C)), 1);
	for c = 1:length(C)
		for v = 1:length(C{c}) 
			VertexNeighbors{C{c}(v)}(end+1) = c;		
		end
	end

end