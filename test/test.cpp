#include <bits/stdc++.h>
#include "test_classes.h"
#include "test_readfile.h"

using namespace std;

void testClasses(){
    cout << "Testing classes" << endl;
    testPrimitive();
    testAtom();
    testBasis();
    testMolecule();
    cout << "Test classes passed" << endl << endl;
}
void testReadFile(){
    cout << "Testing read file" << endl;
    testReadInput();
    //testReadOutput();
    cout << "Test read file passed" << endl << endl;
}

int main(){
    testClasses();
    testReadFile();
}