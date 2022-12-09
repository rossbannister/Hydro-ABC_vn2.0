All alterations to dry ABC-model are commented with ‘Hydro’.

The way reading in namelist is changed to read from standard io in Master_PrepareABC_InitState.f90 and Master_RunNLModel.f90, so that namelist can have any filename.

e.g.  Master_RunNLModel.out < path_to_namelist/name_of_namelist
