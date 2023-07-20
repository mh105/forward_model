Helsinki BEM Framework LCISA solver for MEG/EEG, distribution, v160405

This Helsinki BEM Framework (HBF) package contains a streamlined version of the BEM solver (linear collocation BEM formulated using isolated source approach, LCISA) that was used as a test solver in articles

Stenroos, M., Hunold, A., Haueisen, J., 2014. Comparison of three-shell and simplified volume conductor models in magnetoencephalography. Neuroimage 94, 337–348.

Stenroos, M., Nummenmaa, A., 2016. Incorporating and compensating cerebrospinal fluid in surface-based forward models of magneto- and electroencephalography. PLoS ONE 11(7):e0159595, 2016.

Please cite the 2014 paper, if you use this code for basic MEG modeling (3-shell or 1-shell), and the 2016 paper, if you use this code for EEG or advanced MEG modeling (4-compartment). If (when...) you wish to cite also my method papers that describe, how the computations are done, please see the references at the end of this document.

The solver is aimed for MEG/EEG use, but it works just as well or better for MCG/ECG. There are more accurate BEM approaches (for example Linear Galerkin BEM and Symmetric BEM), but taking into account the level of approximation in the head and thorax models, this solver is adequate for all experimental use, when some attention is paid to the distance of source from nearest boundaries (see the papers).

The example scripts present all necessary file formats. Regarding boundary meshes, it is essential to know that
1) meshes are ordered from inwards to outwards
2) meshes must be closed and
3) triangle orientation is CCW --- please run hbf_CheckTriangleOrientation to check this
4) the resulting models have been verified with meshes of approx 2500 vertices per surface.
5) the recommended distance for dipole sources from the inner skull is half of triangle side length.
6) all physical measures are in SI units --- spatial variables in meters, dipole moments in ampere meters.
7) all sources must be inside the ISA surface --- in the case of three-shell MEG models, inside surface one.

Before using the solver with your own meshes, please study the example case to see, what kind of meshings and parameters to use. If you get surprising results (like essentially worse performance than that presented in the paper), please check the above list and plot your model geometry. If there is something strange, please drop me a line! And, don't worry about the computational cost with these 2500-vertex meshes --- my desktop PC builds the whole model in less than 40 seconds, so there really is no need to use coarser meshes.

The solver is copyrighted. You may use it for academic purposes, but not for anything commercial. For the time being, you are not allowed distribute this package outside our own workgroup, and you are not allowed to put any of these files for public download. These conditions will be relaxed, when I have time for making some kind of www page, documentation and proper distribution package. 

If you wish to include this package as part of your own MEG/EEG toolbox (please do!), please contact me, and we'll sort the necessary things out. I'd be happy to have this solver included in FieldTrip and Brainstorm packages, but I am too busy to do that myself. 

20 Dec 2016 

Matti Stenroos (matti.stenroos@aalto.fi) 
	Department of Neuroscience and Biomedical Engineering
	Aalto University

References

Stenroos, M., Nummenmaa, A., 2016. Incorporating and compensating cerebrospinal fluid in surface-based forward models of magneto- and electroencephalography. PLoS ONE 11(7):e0159595, 2016.

Stenroos, M., Hunold, A., Haueisen, J., 2014. Comparison of three-shell and simplified volume conductor models in magnetoencephalography. Neuroimage 94, 337–348.

Stenroos, M., Sarvas, J., 2012. Bioelectromagnetic forward problem: isolated source approach revis(it)ed. Phys Med Biol 57, 3517–3535.

Stenroos, M., Nenonen, J., 2012. On the accuracy of collocation and Galerkin BEM in the EEG/MEG forward problem. Int J Bioelectromagnetism 14, 29–33.

Stenroos, M., Mäntynen, V., Nenonen, J., 2007. A Matlab library for solving quasi-static volume conduction problems using the boundary element method. Comput Methods Programs Biomed 88, 256–263.

For BEM solution in a general piece-wise homogeneouns model, see

Stenroos M., 2016. Integral equations and boundary-element solution for static potential in a general
 piece-wise homogeneous volume conductor. Phys Med Biol 61:N606–N617.

