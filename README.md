# CM5-calculator
## Description

Allows the calculation of CM5 charges by supplying only the structure and Hirshfeld charges. No need for any additional input.

Based on the charge model 5 (CM5) by Marenich et al. ([10.1021/ct200866d](https://dx.doi.org/10.1021/ct200866d)).
Values inside are copied from their paper or the respective reference they give (CRC Handbook of Chemistry and Physics, 91st ed. (2010–2011) ; 2010; p 9– 49).

Please do contact me with questions or use the GitHub features to contribute.

## Requirements:

* Python version 3 (any should be fine)
* Numpy (any recent version should be fine)
* [ASE](https://wiki.fysik.dtu.dk/ase/index.html)

## Usage:
```python
from cm5-calculator import calc_cm5
from ase import io

molecule = io.read('mol.xyz')
hirshfeld = get_your_charges_somehow()
cm5 = calc_cm5(molecule, hirshfeld)
```
