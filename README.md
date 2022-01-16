<font size = "5">This is a Hartree-Fock program only for s-type orbitals.</font> 

# Setup
The executable files are present so this step can be ignored. </br>
In order to compile the codes, c++ 17 and the two libraries, Eigen and Boost, are required.<br/>

## Install Eigen3
<ol>
<li> Run: <br/> <code>sudo apt install libegin3-dev</code> </li>
<li> Replace EIGENPATH in the makefiles(Makefile, test/Makefile) to the installed path. </li>
</ol>

## Install Boost
<ol>
<li> Run: <br/> <code>sudo apt-get install libboost-all-dev</code> </li>
</ol>

# Testing
<ol>
<li> Enter subdirectory test/:<br/> 
<code>cd test/</code> </li>
<li> Compile: (The executable file is present so this step can be ignored)<br/> 
<code>make</code> </li>
<li> Run: <br/> <code>./test</code> </li>
<li> Check that "All tests passed" is printed on the terminal. </li>
</ol>

# Main Program
<ol>
<li> Return to the main directory:<br/> 
<code>cd ..</code> </li>
<li> Compile: (The executable file is present so this step can be ignored)<br/> 
<code>make</code> </li>
<li> Run: <br/> <code>./run_hf</code> </li>
<li> Enter the input file path, and the log file will be generated at the same directory as the input file.</li>
</ol>

# Diary
## 2022/01/11
<ul>
<li> Prepared environment, installed Eigen3 </li>
<li> Constructed structure of program </li>
<li> Completed classes: Primitive, Basis, Atom </li>
<li> Molecule class 20% done (members) </li>
<li> Unit test for classes </li>
</ul>

## 2022/01/12
<ul>
<li> Update Atom class (added mass) </li>
<li> Update unit test for Atom and Molecule classes </li>
<li> Molecule class 70% done (calculate integrals) </li>
<li> Added integrals.cpp for integral calculation </li>
</ul>

## 2022/01/13
<ul>
<li> Update Atom, Basis, Molecule class (change Atom to member of Basis) </li>
<li> Completed Molecule class </li>
<li> Installed Boost library </li>
<li> Completed read_input.cpp </li>
</ul>

## 2022/01/14
<ul>
<li> Moved integral calculation from Molecule to integral.cpp </li>
<li> Added unit test for read_input </li>
</ul>

## 2022/01/15
<ul>
<li> Bug fix on integral.cpp/_boysFunction() </li>
<li> Added unit test for integrals </li>
<li> Completed HF </li>
</ul>

## 2022/01/16
<ul>
<li> Completed stable version 1.0 </li>
<li> Added unit test for read_output, final energies </li>
</ul>