# Evidence that variation in root anatomy contributes to local adaptation in Mexican native maize
Chloee M. McLaughlin, Meng Li, Melanie Perryman, Adrien Heymans, Hannah Schneider, Jesse R. Lasky and Ruairidh J. H. Sawers

# GRANAR - MECHA
This repository has a copy of the script used to generate the root cross-section anatomies and the estimation of the radial and axial hydraulic properties for the entitled scientific paper.

### Used known model: 
- The Generator of root anatomy in R - **GRANAR** [Heymans et al. 2020](https://doi.org/10.1104/pp.19.00617)
- The model of explicit hydraulic anatomy - **MECHA** [Couvreur et al. 2018](https://doi.org/10.1104/pp.18.01006)

## Input
The GRANAR script used a set of 5 input parameter to get a full cellular network of a root cross-section

| Anatomical parameter | Unit |
|----------------------|------|
| Root radius          | µm   |
| Stele radius         | µm   |
| Aerenchyma area      | µm^2 |
| Metaxylem area       | µm^2 |
| Nbr of metaxylem elm | -    |

The MECHA script used a set of 4 base input parameter

| Sub-cell hydraulic properties   | Unit               |
|---------------------------------|--------------------|
| Cell wall conductivity          | cm^2/hPa/d         |
| Cell membrane base permeability | cm/hPa/d           |
| Aquaporin contribution          | cm/hPa/d           |
| Conductance of plasmodesmata    | cm^3/hPa/d/plasmo. |
| Cell wall thickness             | µm                 |

## Results
...
