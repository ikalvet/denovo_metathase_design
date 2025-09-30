echo "Running an example design job"
python post_rifdock_FastDesign.py --pdb inputs/example_rifdock_output_c0.pdb --nstruct 1 --repack --min_polar 2 --outdir ../results/ --params ligand/HGS.params --ramp_cst_weights 10.0 --no_ala_design --dalphaball ~/DAlphaBall/src/DAlphaBall.gcc
