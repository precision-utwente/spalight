# Validating SPACAR model input

This is a list of requirements for a valid SPACAR model, and actions to take if they are not met.

1. Nodes

    - list of unique, feasible coordinates in space [error]

2. Elements
    
    - unique pairs of existing node numbers (p/q order does not matter) [error]
    - all node numbers should be 'used' at least once [error, but showGeom]
    - element length larger than .. [error]

3. Node properties

    no fixed BC and force/moment/displ/rot in the same direction [error, but showGeom]
    no all-fixed BC and mass/inertia [error, but showGeom]

4. Element properties

    if no flexibility specified within set, no further input needed, but allowed for visualization purposes (color, orien, hide, type, dim) [warn about redundant E,G,rho,nbeams input]
    if an element property set is created, El_Nrs is required [error]
    El_Nrs for all sets combined should be unique and existing
    (not all elements require specification though) [error]
    if orient specified, check if it works. If not specified, check if default works [error]
    if flexibility specified, require E, G, rho, type, dim

    if no flexibility in any set, no need to run spacar (note user though)

5. Releases

6. Simulation properties