# Research data supporting "Mechanical properties of graphene oxide from machine-learning-driven simulations"

This repository contains code and data supporting the following work:

> Zakariya El-Machachi, Bowen Cheng, Volker L. Deringer, _in preparation_

The purpose is to enable readers to access the structural models and simulation code for reproducing the work. 

## Contents

The repository is structured as follows:

* **Code**: The straining code used in this work. This code has been modified from the original production code only in documentation to remove typos and to remove absolute paths. The Young's modulus is calculated as follows: 
    ```
    def get_young_modulus(stress, strain):
        # Fit a linear model to the stress-strain data in elastic region 
        coefficients = np.polyfit(strain[1:101], stress[1:101], 1)
        young_modulus = coefficients[0]
    return young_modulus
    ```
* **Structures**: 10 input structures for GO and 10 for rGO.  

The MACE model used in this work is available at [zakmachachi/GO-MACE-23](https://github.com/zakmachachi/GO-MACE-23).

Due to large file sizes, trajectory data will be shared via Zenodo in due course.

## Figures
The following figures were constructed as follows:
* **Figure 1**: Structural model can be found at [zakmachachi/GO-MACE-23](https://github.com/zakmachachi/GO-MACE-23). Trajectory data will be shared via Zenodo in due course.
* **Figure 2**: Input structures can be found [here](structures). For panel a), b) and c), the structures are in `p1-p2-GO.xyz`. The table below maps the structure configurations for a "batch":

    - **Top row**: OH/O ratio
    - **Left column**: O content (%)
    - **Cells**: Index of each configuration

    | OH/O → <br> O content (%) ↓ | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 |
    |-----------------------------|------|------|------|------|------|
    | 0.10                        | 0    | 1    | 2    | 3    | 4    |
    | 0.20                        | 5    | 6    | 7    | 8    | 9    |
    | 0.30                        | 10   | 11   | 12   | 13   | 14   |
    | 0.40                        | 15   | 16   | 17   | 18   | 19   |
    | 0.50                        | 20   | 21   | 22   | 23   | 24   |

    In total there are 10 batches, meaning 10 structures per config type, and having 25 config types means a total of 250 structures total. Since we are straining in _x_ and in _y_, we plot the average of the combined _x_ and _y_ data meaning we have 20 configurations per line in the plot (shaded with ± 1 standard deviation). 
* **Figure 3**: Input structures can be found [here again](structures). For panel a), the structures are in `p1-p2-GO.xyz`, for panel b) they are in `p1-p2-rGO.xyz`. We are now separating _x_ and _y_ values of $E$, which are slightly different. Each data point is the average $E$ of the 10 config types for _x_ strained (full symbol) and the same 10 strained in _y_ (hollow symbol).  rGO will lose some species under thermal annealing and since these are removed, it results in a lower O% content in some cases. We also see a transformation in functional groups and so this is also in reference to the initial structure.
