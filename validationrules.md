# Validating SPACAR model input

This is a list of requirements for a valid SPACAR model [and actions to take if they are not met].

1. Nodes

2. Elements
    
    - unique pairs of existing node numbers (p/q order does not matter) [error]
    - all node numbers should be 'used' at least once [warning, but showGeom]
    - element length larger than .. [error]

3. Node properties

    - at least six fixes or inputs on the set of nodes [error]
    - no fixed BC and force/moment/displ/rot in the same direction [warning, but showGeom]
    - no all-fixed BC and mass/inertia [warning, but showGeom]

4. Element properties

    - if no flexibility specified within set, no further input needed, but allowed for visualization purposes (color, orien, hide, type, dim) [warn about redundant E,G input]
    - if an element property set is created, El_Nrs is required [warning]
    - El_Nrs for all sets combined should be unique and existing
    (not all elements require specification though) [error]
    - if orient specified, check if it works. If not specified, check if default works [error]
    - if orient specified, check if normal to local x axis [warning]
    - if flexibility specified, require E, G, rho, type, dim [error]

    - if no flexibility in any set, no need to run spacar (note user though) [warning]

5. Releases

    - check if specified releases are already specified as flexible [error]

6. Simulation properties