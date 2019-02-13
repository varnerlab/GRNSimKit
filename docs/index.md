## Gene regulatory network simulation kit (GRNSimKit)
The gene regulatory network simulation kit (GRNSimKit) is a [Julia package](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) to simulate gene regulatory networks using biophysical models with parameters largely derived from primary
literature or parameter databases such a [BioNumbers](https://bionumbers.hms.harvard.edu).

### Installation
GRNSimKit is a [Julia package](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) which can be installed in the ``package mode`` in Julia.
Start of the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/index.html) and enter the ``package mode`` using the ``]`` key (to get back press the ``backspace`` or ``^C`` keys). At the prompt enter:

    (v1.1) pkg> add https://github.com/varnerlab/GRNSimKit.git

This will install GRNSimKit and other all required packages.

### Testing

###

### Code
* [Writing a model file](/model/index.md)
* [How do we load a data dictionary?](/data/index.md)
* [How do we solve a system of GRN equations?](/solve/index.md)
