##
##      This is the startup file for tkmedit when used in conjunction with mne_analyze
##
##      Copyright 2007
##
##      Matti Hamalainen
##      Athinoula A. Martinos Center for Biomedical Imaging
##      Massachusetts General Hospital
##      Charlestown, MA, USA
##
##      No part of this program may be photocopied, reproduced,
##      or translated to another program language without the
##      prior written consent of the author.
## 
##      $Header$
##      $Log$
##      Revision 1.1  2007/04/11 14:55:09  msh
##      tkmedit startup script to be used when launched from mne_analyze
##
##
##
#
#   Pretend that we are interactive
#
set tcl_interactive 1
#
#   Default display configuration
#
SetDisplayConfig 2 2 1
#
#   Show RAS coordinates by default
#
ShowRASCoords 1
#
#   Do not show the program controls by default
#
wm state . withdrawn
#
#   Standard surface loads
#
#
#LoadMainSurface 0 lh.white
#LoadMainSurface 1 rh.white
#
#   Do not show them by default
#
#SetDisplayFlag 4 0
#SetDisplayFlag 5 0
#SetDisplayFlag 6 0


