#include <bits/stdc++.h>
#include "test_classes.h"
#include "../class/class.h"
#include "../src/utils.h"

using namespace std;

void testPrimitive(){
    Primitive p1(3.0, 2.5);
    Primitive p2(2.0, 1.5);

    assert(isEqual(p1.alpha, 3.0) && "Incorrect alpha value for primitive p1");
    assert(isEqual(p1.c, 2.5) && "Incorrect coefficient value for primitive p1");
    assert(isEqual(p2.alpha, 2.0) && "Incorrect alpha value for primitive p2");
    assert(isEqual(p2.c, 1.5) && "Incorrect coefficient value for primitive p2");
    assert(isEqual(getZeta(p1, p2), 5.0) && "Incorrect zeta value between p1 and p2");
    assert(isEqual(getXi(p1, p2), 1.2) && "Incorrect xi value between p1 and p2");
    cout << "Test primitive passed" << endl;
    return;
}
void testBasis(){
    Primitive p1(3.0, 2.5);
    Primitive p2(2.0, 1.5);
    vector<Primitive> p1s{p1, p2};
    Basis b1(p1s);
    p1s.clear();

    assert(isEqual(b1.primitives[0].alpha, 3.0) && "Incorrect alpha value for the first primitive");
    assert(isEqual(b1.primitives[0].c, 2.5) && "Incorrect coefficient value for the first primitive");
    assert(isEqual(b1.primitives[1].alpha, 2.0) && "Incorrect alpha value for the second primitive");
    assert(isEqual(b1.primitives[1].c, 1.5) && "Incorrect coefficient value for the second primitive");
    assert(b1.n_primitive == 2 && "Incorrect n_primtive value for basis b1");
    cout << "Test basis passed" << endl;
    return;
}
void testAtom(){
    Primitive p11(3.0, 2.5);
    Primitive p12(2.0, 1.5);
    vector<Primitive> p1s{p11, p12};
    Basis b1(p1s);
    p1s.clear();

    Primitive p21(2.5, 3.5);
    Primitive p22(0.5, 4.0);
    Primitive p23(3.5, 1.0);
    vector<Primitive> p2s{p21, p22, p23};
    Basis b2(p2s);
    p2s.clear();

    vector<double> coord1{0.3, -0.6, 1.4};
    vector<double> coord2{1.3, 2.0, -1.5};
    vector<Basis> basis_set{b1, b2};
    Atom atom1(coord1, basis_set);
    Atom atom2(coord2, basis_set);
    coord1.clear();
    coord2.clear();
    basis_set.clear();
    assert(isEqual(atom1.basis_set[0].primitives[0].alpha, 3.0) && "Incorrect alpha value for basis_set[0].primitives[0]");
    assert(isEqual(atom1.basis_set[1].primitives[0].c, 3.5) && "Incorrect coefficient value for basis_set[1].primitives[0]");
    assert(atom1.n_basis == 2 && "Incorrect n_basis value for atom1");
    assert(atom2.n_basis == 2 && "Incorrect n_basis value for atom2");
    assert(atom1.basis_set[0].n_primitive == 2 && "Incorrect n_primitive value for atom1.basis_set[0]");
    assert(atom2.basis_set[1].n_primitive == 3 && "Incorrect n_primitive value for atom2.basis_set[1]");
    assert(isEqual(atom1-atom1, 0) && "Incorrect distance between atom1 and atom1");
    assert(isEqual(atom1-atom2, 4.02119385258) && "Incorrect distance between atom1 and atom2");
    assert(isEqual(atom2-atom1, 4.02119385258) && "Incorrect distance between atom2 and atom1");
    cout << "Test atom passed" << endl;
    return;
}
void testMolecule(){
    Primitive p11(3.0, 2.5);
    Primitive p12(2.0, 1.5);
    vector<Primitive> p1s{p11, p12};
    Basis b1(p1s);
    p1s.clear();

    Primitive p21(2.5, 3.5);
    Primitive p22(0.5, 4.0);
    Primitive p23(3.5, 1.0);
    vector<Primitive> p2s{p21, p22, p23};
    Basis b2(p2s);
    p2s.clear();

    vector<double> coord1{0.3, -0.6, 1.4};
    vector<double> coord2{1.3, 2.0, -1.5};
    vector<Basis> basis_set{b1, b2};
    Atom atom1(coord1, basis_set);
    Atom atom2(coord2, basis_set);
    coord1.clear();
    coord2.clear();
    basis_set.clear();
    
    vector<Atom> atoms{atom1, atom2};
    // dummy integrals for mol_explicit
    vector<vector<double>> s_itg;
    vector<vector<double>> t_itg;
    vector<vector<double>> v_itg;
    vector<vector<vector<vector<double>>>> e_itg;
    Molecule mol_implicit(1, atoms);
    Molecule mol_explicit(1, atoms, s_itg, t_itg, v_itg, e_itg);
    atoms.clear();
    assert(isEqual(mol_implicit.atoms[0].basis_set[0].primitives[0].alpha, 3.0) && "Incorrect alpha value for atoms[0].basis_set[0].primitives[0] (mol_implicit)");
    assert(isEqual(mol_implicit.atoms[1].basis_set[1].primitives[0].c, 3.5) && "Incorrect coefficient value for atoms[1].basis_set[1].primitives[0] (mol_implicit)");
    assert(mol_implicit.atoms[0].n_basis == 2 && "Incorrect n_basis value for atom1 (mol_implicit)");
    assert(mol_implicit.atoms[1].n_basis == 2 && "Incorrect n_basis value for atom2 (mol_implicit)");
    assert(mol_implicit.atoms[0].basis_set[0].n_primitive == 2 && "Incorrect n_primitive value for atoms[0].basis_set[0] (mol_implicit)");
    assert(mol_implicit.atoms[1].basis_set[1].n_primitive == 3 && "Incorrect n_primitive value for atoms[1].basis_set[1] (mol_implicit)");
    assert(isEqual(mol_implicit.atoms[0]-mol_implicit.atoms[0], 0) && "Incorrect distance between atom1 and atom1 (mol_implicit)");
    assert(isEqual(mol_implicit.atoms[0]-mol_implicit.atoms[1], 4.02119385258) && "Incorrect distance between atom1 and atom2 (mol_implicit)");
    assert(isEqual(mol_implicit.atoms[1]-mol_implicit.atoms[0], 4.02119385258) && "Incorrect distance between atom2 and atom1 (mol_implicit)");
    assert(mol_implicit.n_atom == 2 && "Incorrect n_basis value for mol_implicit");
    assert(mol_implicit.charge == 1 && "Incorrect charge value for mol_implicit");
    
    assert(isEqual(mol_explicit.atoms[0].basis_set[0].primitives[0].alpha, 3.0) && "Incorrect alpha value for atoms[0].basis_set[0].primitives[0] (mol_explicit)");
    assert(isEqual(mol_explicit.atoms[1].basis_set[1].primitives[0].c, 3.5) && "Incorrect coefficient value for atoms[1].basis_set[1].primitives[0] (mol_explicit)");
    assert(mol_explicit.atoms[0].n_basis == 2 && "Incorrect n_basis value for atom1 (mol_explicit)");
    assert(mol_explicit.atoms[1].n_basis == 2 && "Incorrect n_basis value for atom2 (mol_explicit)");
    assert(mol_explicit.atoms[0].basis_set[0].n_primitive == 2 && "Incorrect n_primitive value for atoms[0].basis_set[0] (mol_explicit)");
    assert(mol_explicit.atoms[1].basis_set[1].n_primitive == 3 && "Incorrect n_primitive value for atoms[1].basis_set[1] (mol_explicit)");
    assert(isEqual(mol_explicit.atoms[0]-mol_explicit.atoms[0], 0) && "Incorrect distance between atom1 and atom1 (mol_explicit)");
    assert(isEqual(mol_explicit.atoms[0]-mol_explicit.atoms[1], 4.02119385258) && "Incorrect distance between atom1 and atom2 (mol_explicit)");
    assert(isEqual(mol_explicit.atoms[1]-mol_explicit.atoms[0], 4.02119385258) && "Incorrect distance between atom2 and atom1 (mol_explicit)");
    assert(mol_explicit.n_atom == 2 && "Incorrect n_basis value for mol_explicit");
    assert(mol_explicit.charge == 1 && "Incorrect charge value for mol_explicit");
    cout << "Test molecule passed" << endl;
    return;
}