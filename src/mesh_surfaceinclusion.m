function [ bounded,in ] = mesh_surfaceinclusion(Faces, Vertices, Inside_Vertices)
%Function used to check whether the surfaces are bounded in that all
%vertices from the inside BEM boundary surface are within the closed
%compartment of the next outer BEM boundary.

in = intriangulation(Vertices,Faces,Inside_Vertices);

if all(in)
    bounded = true;
else
    bounded = false;
end

end

