# CT_SEM_alignment
Programmes to find a thin section in a CT stack using SIFT

Requires python 2 and FIJI

####################################

functions.py 

A header file containing functions to load save and manipulate image stacks in python.

############################

save_sift_correspondenc.ijm

An imageJ macro which compares each slice in a CT stack to the 2D back scattered electron image and save a table of the matching points.

########

plane_fit.py

Contains functions required to fit a plane/paraboloid to the saved from imageJ/

#####



###GUIDE####

1)Make sure the CT stack is the approximate orientation of the thin section, save as an image sequence.

2)Make a directory to save the matching points from fiji too, say "/fiji_out/"

3)Open fiji, and load the macro "save_sift_correspondenc.ijm"

4)Edit the paths to the appropriate names for the current problem

5)Run the script.

6)Make a directory to save the alignment, say "/aligned/", in the directory make another directory called "/aligned_stack/"

7)Open plane_fit.py and edit the paths under __main__ to the suitable directories

8)Run plane_fit.py


###Notes####

Try editing the 'maximal_alignment_error' value in save_sift_correspondenc.ijm to get more or less matches.

Once an idea of where the slice might be, narrow the number of slices you look through in save_sift_correspondenc.ijm.

If SIFT matches cannot be found. Reorientate the stack more accuratley to the direction of the thin section.

#####################
