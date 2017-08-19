# Syntax

## General
- A simulation is performed with the spacarlight() command.
- The spacarlight() command requires at least five input arguments (a sixth is optional). If fewer are supplied, warnings and a visualization can appear to help in completing the input.
- The results of a successful simulation are returned in a structure (the first and only output argument)

## Input
1. Node positions nodes

2. Element connectivity elements

3. Node properties nprops.

    fix
    fix_pos 
    fix_orien
    
    displ_x
    displ_y
    displ_z
    rot_x
    rot_y
    rot_z

    force
    moment

    mass
    mominertia

    fix_x
    fix_y
    fix_z

    force_initial (force_add)
    moment_initial (moment_add)

    displ_initial_x (displ_x_add)
    displ_initial_y (displ_y_add)
    displ_initial_z (displ_z_add)
    displ_initial_rx (rot_x_add)
    displ_initial_ry (rot_y_add)
    displ_initial_rz (rot_z_add)
    
4. Element properties eprops.

    elems
    emod
    smod
    dens
    type
    dim
    orien
    nbeams
    flex
    color
    hide

5. Releases rls.

    def

6. Optional input opt.
    
    filename
    gravity
    silent
    calcbuck
    showinputonly