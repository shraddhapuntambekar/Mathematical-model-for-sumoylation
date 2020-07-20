# <b>Mathematical model for sumoylation</b>
<b>Matlab codes to simulate mechanism-based, ODE model for a multi-step, multi-enzymatic post-translational modification in which a small ubiquitin-like modifier (SUMO) is attached to the target. The model includes enzyme mechanism details such as autosumoylation of E2 and also the bifunctional nature of SENP.</b>

This is the first mathematical model for the sumoylation systems, the results of which are published in Puntambekar, S. S., Nyayanit, D., Saxena, P., & Gadgil, C. J. (2016). Identification of Unintuitive Features of Sumoylation through Mathematical Modeling. Journal of Biological Chemistry, 291(18), 9458-9468 (https://www.jbc.org/content/291/18/9458.full).
Codes included in this project are to stimulate the system's concentration dynamics, to study steady state analysis (numerically and analytically) and to examine the effect of varying system's parameters on its steady state. 
This model can be used to generate testable hypotheses for existence and mechanisms of unintuitive behaviours.

<b>(1) Folder: timecourse_steadystate/ </b>    
The code allows to simulate a large range of experimental conditions, each being one specific combination of the following multiple possible scenarios: 
- Closed system: A phenomena occuring at short time scales; where there is no protein synthesis and degradation can be simulated. This case represents an in-vitro system where short term responses are monitored.
- Open system : Long term events, with zero order formation of enzymes and targets and their first order degradation, can be simulated. These represents an in-vivo experiment where long term response is monitored.
- E2 as target: One of the interesting features of sumoylation system is modification of its own enzyme (E2) along with the targets. Sumoylation of E2 might alter its structure and hence have a different binding affinity towards a target then the unmodified E2. The code allows to simulate a system with or without E2's modification.
- System having 2 targets: Code allows to simulate a system with two targets, each having different binding affinities towards (un)sumoylated E2.
- time course analysis: Code allows to simulate and track the changes in concentrations of the system's components over time. 
- steady state analysis : Code allows to simulate and predict steady state concentrations of the system. Steady state is numerically calculated using matlab's function 'fsolve'.

This folder contains 10 matlab files:
- main.m: The main file to run the code and to simulate large range of possible experimental conditions of the sumoylation system
- getnames.m: the script where the names of the system's components (enzymes and targets) are defined
- getinitial.m: script to define initial concentrations of the system's components
- getinitialCLOSED.m: script to define initial concentrations of the proteins and enzymes; for a closed system
- getparams_real.m: script to assign values for all the model parameters. The values assigned currently are approximated from enzyme kinetics experiments as tabulated in are paper and also documented in the file 'parameter_values_usedinmodel_with_refs.xlsx'.
- getrate.m: script where the reaction rates for each step of the sumo cycle are defined. 
- getrate_forfsolve.m: Function with algebraic equations of an open system obtained on equating all the derivatived to zero.
- getrate_forfsolve_closed.m: Function with algebraic equations of a closed system obtained on equating all the derivatived to zero. 

<b>(2) Folder: analytical_steadystate/ </b>  
This folder contains a mathematica script for analytical steady state analysis of the sumoylation system.

<b>(3) parameter_variation/ </b>   
Code to vary the system's parameters and examine their effect on the system's steady state. The following 3 systems parameters are varied: 
- parameter for SENP formation  (param.k(25,1)) - As SENP has two seemingly opposite functions - of both sumoylation and desumoylation; it is intuitively difficult and hence interesting to predict the effect of change in SENP levels on target's sumoylation.
- binding affinity of modified and unmodified E2 to the target (param.k(14,1),param.k(14,2)) -It is experimentally observed that sumoylation of E2 can either increase, decrease or have no effect on target sumoylation. We simulate cases where we vary the target bindin parameter of sumoyalted E2 (param.k(14,2)) and that of unsumoylated E2 param.k(14,1).
- parameters for degradation of sumoE1 and sumoE20 (param.k(27,1), param.k(28,1))  - On analytically solving a case of simple open system with sumoylation of only one target and no autosumoylation of E2, we observed that target modification is robust to properties of modifying enzyme E2. Hence its interesting to study the effect of varying the degradation parameters of the intermediate complexes.

This folder contains 7 matlab scripts with 'main_paravar_e2notTar_1tar_open.m' being the main file to  run the scripts and simulate the system.

<b>(4) mTar_nSumo/ </b>   
This folder contains scripts to simulate a simple sumo system in which user can define multiple number of targets each having different levels of polysumoylation. This modeleled system is slightly different than the one whose results we have discussed in the paper mentioned above. This system considers sumoylation of another of the system's enzyme E1 along with E2. 

<b>(5) ModelAssumptions.txt </b>  
This text file contains a list of assumptions that were made while constructing the mathematical model for sumoylation system.

<b> (6) SumoSytem.jpg </b>  
Pictorial representation of the modeled system.

<b>(7) parameter_values_usedinmodel_with_refs.xlsx </b>   
The file contains values for the reaction rate parameters used in the simulations each citing a reference from literature.







  
 
