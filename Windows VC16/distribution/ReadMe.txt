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

 - a tutorial on using the MESMER GUI (a web program) to construct simple data files.

To check that the program is working properly, open a DOS window in the Mesmer folder,
move to its MesmerQA sub-folder, and type QA.bat.
This runs several example files and checks that the output is as expected. 

To edit and recompile the C++ source code for Mesmer, download the project from
SourceForge via SVN.

Mesmer on SourceForge: https://sourceforge.net/projects/mesmer/
   
----------------------------------------------------------------------------------------------------
* Struan H. Robertson Univ. Leeds    struanhrobertson@gmail.com
* Paul Seakins        Univ. Leeds    P.W.Seakins@leeds.ac.uk
* David Glowacki      Univ. Bristol  David.R.Glowacki@bristol.ac.uk
* Robin Shannon       Univ. Leeds    shannon.rj@gmail.com
* Mark Blitz          Univ. Leeds    M.A.Blitzg@leeds.ac.uk

