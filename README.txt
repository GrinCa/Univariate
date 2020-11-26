*****************************************************************************************
#################################       USEFUL        ###################################
*****************************************************************************************
In case you only want the source code, please go to https://github.com/GrinCa/Univariate 
and clone the repository with git:

git clone https://github.com/GrinCa/Univariate

Do not hesitate to contact me (Julien Ecobichon) for any questions, I would be happy to 
help you (Matlab/Gmsh/FreeFem++):
-email: julien9844.je@gmail.com

It can be useful to look at my report: "Efficient modelling of acoustic problems". You 
can either ask me (email) or ask M. Rumpler to get it.

In case your work is to improve this code, I recommend you to use software for versionning
like git. 

*****************************************************************************************
#################################      IMPORTANT        #################################
*****************************************************************************************

this code requires the use of Gmsh and FreeFem++.

This program is not a standalone and therefore is likely to require specific 
modifications when creating a new study, which is by the way the most inte-
resting part!

I strongly advise to download the source code on github. Then you know that you have the
latest version which should be fine.

It's really likely that the code doesn't compile because one file is missing, in this case
do not hesitate to ask me.

*****************************************************************************************
#################################      INFORMATION      #################################
*****************************************************************************************

This code has been designed for univariate studies for the following system:

A(f)U(f) = F(f)     EQUATION 1

where A(f) and F(f) are calculated from the finite element method. Note that
A must be written as A = sum( gi(f)*Ai ) where gi(f) are scalar functions and 
Ai DOESN'T DEPEND of f. Be careful in case Ai depends of f ( it's the case if we 
use PML). Some work have been made in order to counter this specific case without 
any guarantee!

This program implements the Univariate algorithm developped by Slone & al.(see 
https://www.researchgate.net/publication/3017281_Well-conditioned_asymptotic_
waveform_evaluation_for_finite_elements). 

Source code of the Slone algorithm is located under the WCAWE folder (WCAWE_basis.m).


The numerical parameters of the study must be entered in the main_multi.m file. 
There are 2 structs for entering the parameters, param and flag. flag is important 
for post-treatment.

Note that this program can run any kind of studies as long as it respects the pattern
of the EQUATION 1. In case you need to study coupled physical systems, you'll need to
look more in details into the Matrices folder (in particular get_matrices.m) where 
you'll need to create a pattern.m file.

An interval strategy has been developped in order to optimize the method, and require 
to pay attention to some details. Those details are developped in my report, but 
information can be given:
- the sub-basis created by the WCAWE must contain 20 vectors (better result / time)
- the length of the sub-intervals is about 20 to 40 Hz with 20 vectors. The risk of 
extending the length is to obtain false result wthout even knowing it.

Some routines have also been developped to post process the results. It's recom-
mended to pay attention to the function convertGEO2VTK inside the folder DataMap. 
This function take the solution vector U named SOL and other arguments to create 
a vtk file.

In case the user need to implement its own model, with A = sum( gi(f)*Ai ), he will 
need to change the values inside coeff_LHS and coeff_RHS in the main_multi.m. Again 
it is recommended to give a look at the file create_cross_derivatives in the folder 
Derivatives, which has been detailed in particular the recursion aspects.


*****************************************************************************************
#################################      SOME ADVISES      ################################
*****************************************************************************************

If you need to create a new model (Gmsh-FreeFem++), I strongly advise you to compile and
debug your .edp file (FreeFem++). Once you have a good model, you can perform sweeps and 
pay attention to the Matlab interface.







