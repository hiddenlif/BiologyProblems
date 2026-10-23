# Executed findings: Biology problems

Three original synthetic problems grounded in cell-therapy assay mathematics. These are educational reference solutions, not product specifications.

## Observed results

At receptor = 100 nM, total ligand = 10 nM and Kd = 1 nM, bound complex is 9.89024 nM and occupancy is 9.89024%, not the 90.909% obtained by substituting total ligand for free ligand. The bead example yields 15 cells/µL with conditional 95% interval [12.2763, 18.1931]. With prevalence 0.0001 and sensitivity 0.95, the no-background 95% detection calculation requires 31,533 cells. At 30,000 cells and a per-cell blank false-positive rate of 0.0001, an any-positive rule has a 95.022% blank false-alarm probability. Requiring at least seven positive events reduces the blank rate to 3.3501%, with power 36.9569% under the specified alternative.

## Scope of evidence

Binding assumes one-to-one equilibrium and consistent concentration units, without avidity or kinetic effects. The bead interval assumes independent Poisson sampling and known bead addition/sample volume; it excludes pipetting, mixing, recovery and viability uncertainty. Rare-cell calculations assume independent Bernoulli classifications with fixed sensitivity/background. These calculations do not establish an empirical limit of detection, a prevalence estimate or a product release criterion.

The accompanying validation receipt identifies the analysis and dependency hashes and reports pytest and notebook execution. Re-run after changing code, sources or dependencies.
