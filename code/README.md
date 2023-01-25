# code for modsim paper "Modelling aspects of the effect of stigma on the prevalence of anxiety and/or depression"


tldr; run `stigma.m`

## important files

`inputs.csv` are all the parameter value, with rows 
	1. variable name, 
	2. min val, 
	3. max val, 
	4. expected val, 
	5. flag for which distribution to use (0 for triangular, 1 for uniform, 2 to remove)


`stigma.m` is the master file which calls any other required functions or files, and produces the plots for the paper 

`lhsonly.m` generates the latin hypercube sampling across all parameters with flags "0" or "1" (technically anything not "2") in the flags from `inputs.csv`. This results in the file `trilhspoints.txt` for the varying params and `fixedParams.txt` for the constant ones (needed by ODE solver in `runlhs` function in `stigma.m`).

`prcconly.m` calculates the multivariate sensitivity using partial rank correlation coefficient, which one can see in the tornado plot. This uses an output of interest of prevalence (1-Un-Us).

`tornadoplot.m` is code to generate a tornado plot. 


## outputs generated

`trilhspoints.txt` are all the parameters being varied for the PRCC sensitivity analysis, from `lhsonly.m`

`fixedParams.txt` are from `lhsonly.m` and are any constant parameters needed by the ODE solver in `runlhs` but not being varied for the sensitivity analysis.

`lhsoutput.mat` is generated in `runlhs` function sitting in `stigma.m` and are saves of the 9th ODE in the system, currently the cumulative number of new/spontaneous "infections" [Apr 2022]. That is, the "outputs of interest" for the PRCC analysis.

`TornadoPlot` .eps and .fig are both outputs from `tornadoplot.m`, depicting the PRCC analysis.

A bunch of other `*.mat` and `*.fig` files are generated locally to reduce time rerunning outputs. 

