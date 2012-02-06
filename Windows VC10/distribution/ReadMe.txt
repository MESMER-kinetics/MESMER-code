           MESMER (Master Equation Solver for Multi Energy-well Reactions) 

models the interaction between collisional energy transfer and chemical reaction
for dissociation, isomerization and association processes.

It is capable of solving the energy grained master equation (EGME) for a unimolecular
system composed of an arbitrary number of wells, transition states, sinks, and reactants.

Mesmer is a command line program and should be run by first opening a DOS window,
normally with the current directory as the folder containing the input file. Then type:

  mesmer infile.xml
  
  
A Start Menu Programs item contains links to:

 - The Mesmer folder. Subfolders contain examples and it is recommended that you base
   your own initial datafiles in parallel subfolders.

 - MESMER manual.pdf

 - two tutorials, on using Gaussian output and add library molecules.

To check that the program is working properly, open a DOS window in the Mesmer folder,
move to its MesmerQA sub-folder, and type QA.bat.
This runs several example files and checks that the output is as expected. 

To edit and recompile the C++ source code for Mesmer, download the project from
SourceForge via SVN.

Mesmer on SourceForge: https://sourceforge.net/projects/mesmer/
   
----------------------------------------------------------------------------------------------------
* Struan H. Robertson                   struanhrobertson@gmail.com
* Chi-Hsiu Liang                        alvyn.liang@gmail.com
* David Glowacki         Univ. Bristol  David.R.Glowacki@bristol.ac.uk
* Chris Morley                          c.morley@gaseq.co.uk
* Michael Pilling        Univ. Leeds    M.J.Pilling@leeds.ac.uk

