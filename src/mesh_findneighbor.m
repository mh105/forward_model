function [vertex_index, face_index] = mesh_findneighbor(f, v_list, n)
% Function used to find the faces and vertices starting from a center
% vertex after n iterations. 

assert(n >= 1, 'Too few iterations specified. Cannot proceed.')

for ii = 1:n
    face_index = any(ismember(f, v_list), 2);
    involved_vertices = f(face_index, :);
    v_list = unique(involved_vertices);
end

vertex_index = v_list;
face_index = find(face_index==1);

end

