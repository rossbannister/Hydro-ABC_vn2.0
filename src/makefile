#================================================================================
# Compile and link options for ubuntu
#================================================================================
#CompileOpts=-fcheck=bounds
CompileOpts=-O3 -std=legacy
LinkOpts=-L/usr/local/lib -lfftpack5.1d -L/usr/lib/x86_64-linux-gnu  -lnetcdff
LinkOptsEig=-L/usr/local/lib -lfftpack5.1d -llapack -ltmglib -lrefblas -L/usr/lib/x86_64-linux-gnu  -lnetcdff

#================================================================================
# Compile and link options for Met Dept RACC cluster
#================================================================================
#CompileOpts=-O3
#LinkOpts=-L/usr/local/lib -L/home/users/sws98rnb/Lib/lib -lnetcdff -lfftpack5.1d
#LinkOptsEig=-L/usr/local/lib -L/home/users/sws98rnb/Lib/lib -lnetcdff -lfftpack5.1d -llapack -ltmglib -lrefblas



# Linking stage
Master_PrepareABC_InitState.out: Master_PrepareABC_InitState.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Write_state_2d.o BoundaryMod.o Boundaries.o Functions.o Energy.o Calc_geost.o Calc_hydro.o Effective_buoyancy.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Deallocate_dims.o Deallocate_um_data.o Deallocate_model_vars.o Moist_Param.o  # JZ
	f95 -o $@ $^ $(LinkOpts)

Master_RunNLModel.out: Master_RunNLModel.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_Averages.o Moist_Param.o  # JZ
	f95 -o $@ $^ $(LinkOpts)

Master_Linear_Analysis.out: Master_Linear_Analysis.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o
	f95 -o $@ $^ $(LinkOptsEig)

Util_WaveSpectra.out: Util_WaveSpectra.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriverWave.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_Averages.o Write_one_field.o Moist_Param.o
	f95 -o $@ $^ $(LinkOptsEig)


# Compiling stage
Master_PrepareABC_InitState.o: Master_PrepareABC_InitState.f90 DefConsTypes.o
	f95 -c Master_PrepareABC_InitState.f90 $(CompileOpts)

DefConsTypes.o: DefConsTypes.f90
	f95 -c $< $(CompileOpts)

SetOptions.o: SetOptions.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_um_data.o: Initialise_um_data.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_dims.o: Initialise_dims.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_model_vars.o: Initialise_model_vars.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Read_um_data_2d.o: Read_um_data_2d.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Set_grid.o: Set_grid.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Set_ht_dep_cons.o: Set_ht_dep_cons.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Process_um_data.o: Process_um_data.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Write_state_2d.o: Write_state_2d.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Read_state_2d.o: Read_state_2d.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

BoundaryMod.o: BoundaryMod.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Boundaries.o: Boundaries.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Functions.o: Functions.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Energy.o: Energy.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Calc_geost.o: Calc_geost.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Calc_hydro.o: Calc_hydro.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Master_RunNLModel.o: Master_RunNLModel.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

ABC_NL_ModelDriver.o: ABC_NL_ModelDriver.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Calc_vert_mom_source.o: Calc_vert_mom_source.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Calc_horiz_div.o: Calc_horiz_div.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Calc_horiz_vort.o: Calc_horiz_vort.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Effective_buoyancy.o: Effective_buoyancy.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

ABC_NL_model.o: ABC_NL_model.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_Averages.o: Initialise_Averages.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Lscales_from_fft.o: Lscales_from_fft.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Magnitude_rms.o: Magnitude_rms.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_ObsSpec.o: Initialise_ObsSpec.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Master_Linear_Analysis.o: Master_Linear_Analysis.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Write_one_field.o: Write_one_field.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Initialise_CVT.o: Initialise_CVT.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Deallocate_dims.o: Deallocate_dims.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Deallocate_Averages.o: Deallocate_Averages.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Deallocate_um_data.o: Deallocate_um_data.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Deallocate_model_vars.o: Deallocate_model_vars.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Util_WaveSpectra.o: Util_WaveSpectra.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

ABC_NL_ModelDriverWave.o: ABC_NL_ModelDriverWave.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

RMS_diff.o: RMS_diff.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

Moist_Param.o: Moist_Param.f90 DefConsTypes.o
	f95 -c $< $(CompileOpts)

all: Master_PrepareABC_InitState.out Master_RunNLModel.out Master_Linear_Analysis.out Util_WaveSpectra.out

clean:
	rm -rf *.o *.out *.mod
