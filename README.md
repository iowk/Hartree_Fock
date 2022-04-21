<font size = "5">This is a Hartree-Fock program only for s-type orbitals.</font> 

# Setup
In order to compile the codes, a c++ 17 compiler and two libraries, Eigen and Boost, are required.<br/>

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
