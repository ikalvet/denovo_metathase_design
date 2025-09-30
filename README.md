This repository presents the computational workflow used to design de novo artificial metathases using Rosetta, as presented in the manuscript:

**"De Novo Design and Evolution of an Artificial Metathase for Cytoplasmic Olefin Metathesis"<br>**
https://www.researchsquare.com/article/rs-5849532/v1
<br>
<br>

**0) Download and install Rosetta with RifGen and RifDock.**

For getting started with installing rifgen and rifdock, download it from here and follow the instructions:<br>
https://github.com/rifdock/rifdock

It requires the use a specific Rosetta build which can be downloaded from here:<br>
https://files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/rosetta_versions/rosetta_source_for_rifdock_hdf5.tar.gz

Compile and install the Rosetta and rifdock software following their corresponding instructions.

Set environment variables: adjust the $WDIR path to the actual path of the downloaded repository, and set `$RIF_PATH` to the actual install directory of your copy of RifDock.<br>
```
export WDIR=/path/to/this/repository
export RIF_PATH=/path/to/rifdock
```

Set up a python environment that contains `numpy` and `pyrosetta`.
For example, this conda environment definition will be more than enough:<br>
https://github.com/ikalvet/heme_binder_diffusion/blob/main/envs/diffusion.yml<br>
```
conda env create -f diffusion.yml
conda activate diffusion
```
<br>


**1) Generate a rotamer interaction field (RIF) around each of the rotamers of the ligand using RifGen.**


Then run `rifgen` on a single ligand conformer using this command:<br>
`$RIF_PATH/rifdock/latest/rifgen @$WDIR/rifgen.flag -extra_res_fa $WDIR/ligand/rotlib/HGS_RR_c0/HGS.params -rifgen:target $WDIR/ligand/rotlib/HGS_RR_c0/HGS.pdb`

And for all ligand conformers generate `rifgen` commands like this:<br>
`for f in $WDIR/ligand/rotlib/HGS_RR_c* ; do nam=$(basename ${f}) ; echo "$RIF_PATH/rifdock/latest/rifgen @$WDIR/rifgen.flag -extra_res_fa ${f}/HGS.params -rifgen:target ${f}/HGS.pdb -rifgen:outdir output_${nam} > output_${nam}.log" >> commands_rifgen ; done`<br>
And then run these commands separately.

**2) Place the interacting rotamers and the ligand into de novo designed protein scaffolds using RifDock.**

First, update the `code/rifdock.flag` file in the section marked with `############## YOUR RifGen OUTPUT ###############` with the corresponding content from the output of your `RifGen` run STDOUT (`output_{nam}.log`). These are flags that point `rifdock` to a specific `rifgen` output dataset.

Input scaffold can be obtained from the PDB, for example under the id 4YXX.

Then run rifdock on an input scaffold `$WDIR/inputs/scaffold_4yxx.pdb`:<br>
`$RIF_PATH/rifdock/latest/rif_dock_test @$WDIR/rifdock.flag -rif_dock:scaffolds scaffold.pdb -rif_dock:scaffold_res $WDIR/inputs/positions.pos -rif_dock:outdir {outdir}`

RifDock needs to be run separately for each combination of input scaffold and ligand conformer RIF dataset.

**3) Design the amino acid sequence around the ligand using Rosetta FastDesign.**

For each RifDock output, use a corresponding command to run a PyRosetta script to design the amino acid sequence around the docked ligand.<br>
`python $WDIR/post_rifdock_FastDesign.py --pdb $WDIR/inputs/example_rifdock_output_c0.pdb --nstruct 5 --repack --min_polar 2 --outdir outputs/ --params $WDIR/ligand/HGS.params --ramp_cst_weights 10.0 0.0 --no_ala_design`

If you get errors related to DALphaBall then please download and compile your own version from https://github.com/outpace-bio/DAlphaBall.git.

The design script produces a scorefile in the output directory. The final designs can be selected by filtering the scorefile using the metrics highlighted in section "Computational design of Ru1·dnTRP" in the Supplementary Materials of the manuscript.

