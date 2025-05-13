# Effect_of_Nanoparticle_half_life_on_delivery_efficiency

Repository for all code regarding the paper "The effect of nanoparticle half-life on tumour delivery efficiency".

# Script for nanoparticle half-life calculation
The file "Constrained_model_documented.py" does the following:
- calculates half-life from experimental data using the two-compartment model.
- the model constrains the y0 (%ID/g) = 100, based on the assumption at t = 0 mins, 100% of the injected dose is in the blood.

All python processing was performed in an Anaconda environment.


# Endothelial cell uptake simulation
The file "Endothelial_cell_uptake_simulation" does the following:
- Visualizes the simulation of nanoparticles into endothelial cells.
- The simulation constrains the entry sites within the endothelial cells.

The simulation was executed in matlab.
