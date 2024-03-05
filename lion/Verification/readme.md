# Verifications

These verifications are used to test that Lion is correctly calculating the dynamics of ions in ion traps.
The physics behind each verification is described in the accompanying paper, along with an explanation of what each verification does.

## List of verifications

| Verification             | Notes |
|--------------------------|-------|
| SecularMotionFrequencies | Calculates secular frequencies of oscillation for a single ion in an ion trap, over a range of different parameters. `SecularFrequencies_fullRF.m` performs the simulations using the full electric field, while `SecularFrequencies_pseudoPot.m` performs the simulations using the pseudopotential approximation for the trap. `PlotPaperGraph.m` plots the graph used in the paper. |
| EquilibriumSeparation    | Uses a simulation to determine the minimum energy positions of ions in a linear chain, and compares these values to published data in the literature. `EquilibriumSeparation.m` runs the simulations, `PlotPaperGraph.m` plots the graph used in the paper. |
| NormalModes_Linear       | Uses a simulation to determine the normal mode spectrum of ions in a linear chain. `NormalModes_Linear.m` runs the simulations, `PlotPaperGraph.m` plots the graph used in the paper. |
| CoolingMechanisms/LangevinBath | Simulates groups of non-interacting atoms that are coupled to a Langevin bath, and shows the atoms thermalize to the correct temperature over the expected time. `LangevinBath.m` performs the calculation, `PlotPaperGraph.m` plots the graph used in the paper. |
| CoolingMechanisms/LaserCool | Simulates a number of atoms with laser cooling applied along one direction, to show that the velocity parallel to the beam is damped. `LaserCool.m` runs the simulations. |


