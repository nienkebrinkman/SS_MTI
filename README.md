# SS_MTI
Single-Station Moment Tensor Inversion package.

## Usage
The source inversions are currently developed to be applied to data from the InSight mission focussed on marsquakes that occured on Sol 173, 183 and 235.
The inversion method is parallelised.
Currently, only the grid-search and direct inversion are available. Metropolis hastings method will follow soon.

Instaseis is used as a forward modeller. A reflectivity code as forward modeller will follow soon.

--Example scripts--
The /Scripts/ folder consists of python scripts that can be used to perform an inversion:
"SourceInversion_script.py" 
This is the most general script that defines a step-by-step guidance to perform a source inversion.
"Interactive_SI_script copy.py"
This script is used to perform the source inversions of the marsquakes and includes plotting routines.
The plotting routines can be found below """Post-Processing"""
"Synthetic_SI_Script.py"
This script is used to perform the source inversions of synthetically created events and includes plotting routines.
The structure is essentially the same as decribed in "Interactive_SI_script copy.py"

## License

Distributed under the MIT License. See `LICENSE` for more information.


## Contact

Nienke Brinkman - nienke.brinkman@erdw.ethz.ch

Project Link: [https://github.com/nienkebrinkman/SS_MTI](https://github.com/nienkebrinkman/SS_MTI)
