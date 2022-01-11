#include <bits/stdc++.h>
#include "atom.cpp"

using namespace std;
class Molecule
{
    /* class Molecule:
    Molecule object with Atoms
    

    Members:
        Variable         Data type                                  Description
        n_atom           int                                        number of atoms
        atoms            vector<Atom>(n_atom)                       atoms
        charge           int                                        charge
        n_basis          int                                        number of basis functions
        s_itg            vector<vector<double>>(n_basis, n_basis)   overlap integrals
        t_itg            vector<vector<double>>(n_basis, n_basis)   kinetic integrals
        v_itg            vector<vector<double>>(n_basis, n_basis)   nucleus-electron coulomb integrals
        e_itg            vector<vector<vector<vector<double>>>>(n_basis, n_basis, n_basis, n_basis)   two-electron integrals
        s_itg_4idx       vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   overlap integrals (4-indexed)
        k_factor         vector<vector<vector<vector<double>>>>(n_basis, n_prim, n_basis, n_prim)   prefactors k for two-electron integrals
        is_s_itg_4idx_calculated    bool                            Record if s_itg_4idx calculated            
    */  
    public:
    int n_atom;
    vector<Atom> atoms;
    int n_basis;
    vector<vector<double>> s_itg;
    vector<vector<double>> t_itg;
    vector<vector<double>> v_itg;    
    vector<vector<vector<vector<double>>>> e_itg;
    vector<vector<vector<vector<double>>>> s_itg_4idx;
    vector<vector<vector<vector<double>>>> k_factor;
    bool is_s_itg_4idx_calculated = false;

    Molecule(int charge,
    vector<Atom> const& atoms_inp, 
    vector<vector<double>> const& s_itg_inp, 
    vector<vector<double>> const& t_itg_inp, 
    vector<vector<double>> const& v_itg_inp, 
    vector<vector<vector<vector<double>>>> const& e_itg_inp){
        /*
        Read integral from input file explicitly
        Input:
            Variable         Data type                                  Description
            charge_inp       int                                        charge
            atoms_inp        vector<Atom>(n_atom)                       atoms
            s_itg_inp        vector<vector<double>>                     overlap integrals (optional)
            t_itg_inp        vector<vector<double>>                     kinetic integrals (optional)
            v_itg_inp        vector<vector<double>>                     nucleus-electron columb integrals (optional)
            e_itg_inp        vector<vector<vector<vector<double>>>>     two-electron integrals (optional)
        */ 
        _initialize(atoms_inp, charge);
        s_itg = s_itg_inp;
        t_itg = t_itg_inp;
        v_itg = v_itg_inp;
        e_itg = e_itg_inp;
    }    
    Molecule(int charge, vector<Atom> const& atoms_inp){
        /*
        Compute integrals implicitly
        Input:
            Variable         Data type                                  Description
            charge_inp       int                                        charge
            atoms_inp        vector<Atom>(n_atom)                       atoms
        */ 
        _initialize(atoms_inp, charge);
        _calcSItg4idx();
        _calcSItg();
        _calcTItg();
        _calcVItg();
        _calcEItg();
    }
    ~Molecule(){
        atoms.clear();
        s_itg.clear();
        t_itg.clear();
        v_itg.clear();
        e_itg.clear();
        s_itg_4idx.clear();
        k_factor.clear();
    }
    double getSItg(int idx1, int idx2){
        /*
        get overlap integral value between wavefunction idx1 and idx2
        Input:
            Variable         Data type                                  Description
            idx1             int                                        the first index
            idx2             int                                        the second index
        Return:
            Variable         Data type                                  Description
            val              double                                     integral value
        */
        if(idx2 > idx1) swap(idx1, idx2)
        double val = s_itg[idx1][idx2];
        return val;
    }
    double getTItg(int idx1, int idx2){
        /*
        get kinetic energy integral value between wavefunction idx1 and idx2
        Input:
            Variable         Data type                                  Description
            idx1             int                                        the first index
            idx2             int                                        the second index
        Return:
            Variable         Data type                                  Description
            val              double                                     integral value
        */
        if(idx2 > idx1) swap(idx1, idx2)
        double val = s_itg[idx1][idx2];
        return val;
    }
    double getVItg(int idx1, int idx2){
        /*
        get nucleus-electron coulomb integral value between wavefunction idx1 and idx2
        Input:
            Variable         Data type                                  Description
            idx1             int                                        the first index
            idx2             int                                        the second index
        Return:
            Variable         Data type                                  Description
            val              double                                     integral value
        */
        if(idx2 > idx1) swap(idx1, idx2)
        double val = s_itg[idx1][idx2];
        return val;
    }    
    double getEItg(int idx1, int idx2, int idx3, int idx4){
        /*
        get two-electron integral value between wavefunction idx1 and idx2
        Input:
            Variable         Data type                                  Description            
            idx1             int                                        the first index
            idx2             int                                        the second index
            idx3             int                                        the third index
            idx4             int                                        the forth index
        Return:
            Variable         Data type                                  Description
            val              double                                     integral value
        */
        if(idx2 > idx1) swap(idx1, idx2)
        if(idx4 > idx3) swap(idx3, idx4)
        double val = itg[idx1][idx2][idx3][idx4];
        return val;
    }

    private:
    void _initialize(vector<Atom> const& atoms_inp, int charge){
        // common initialization for explicit and implicit integral
        atoms = atoms_inp;
        n_atom = atoms.size();        
        assert(n_atom > 0 && "Number of atoms should be positive");
        charge = charge_inp;
        _calcNBasis();
        return;
    }
    void _calcNBasis(){
        // calculate total number of basis functions
        n_basis = 0;
        for(const auto& atom: atoms){
            n_basis += atom.n_basis;
        }
        return;
    }
    void _calcSItg4idx(){
        // calculate overlap integrals s_itg and 4-indexed overlap integral s_itg_4idx
        is_s_itg_4idx_calculated = true; // set s_itg_4idx as calculated
        for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1){
            
            atom1 = atoms[atom1_idx];
            for(const auto& basis1: atom1.basis_set){
                for(int atom2_idx = 0 ; atom2_idx <= atom1_idx ; ++atom2){
                    atom2 = atoms[atom2_idx];
                    for(const auto& basis2: atom2.basis_set){

                    }
                }
            }
        }
        return;
    }
    void _checkSItg4idxCalculated(){
        assert(is_s_itg_4idx_calculated && "s_itg_4idx is not calculated, should call _calc_s_itg_4idx() first");
        return;
    }
    void _calcSItg(){
        // calculate overlap integrals s_itg
        _checkSItg4idxCalculated()
        return;
    }
    void _calc_TItg(){
        _checkSItg4idxCalculated()
        return;
    }
    void _calc_VItg(){
        _checkSItg4idxCalculated()
        return;
    }
    void _calcKFactor(){
        return;
    }
    void _calcEItg(){
        _calcKFactor();
        return;
    }
};
Molecule::Molecule(int charge,
vector<Atom> const& atoms_inp, 
vector<vector<double>> const& s_itg_inp, 
vector<vector<double>> const& t_itg_inp, 
vector<vector<double>> const& v_itg_inp, 
vector<vector<vector<vector<double>>>> const& e_itg_inp){
    /* Constructor
    Read integral from input file explicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
        s_itg_inp        vector<vector<double>>                     overlap integrals (optional)
        t_itg_inp        vector<vector<double>>                     kinetic integrals (optional)
        v_itg_inp        vector<vector<double>>                     nucleus-electron columb integrals (optional)
        e_itg_inp        vector<vector<vector<vector<double>>>>     two-electron integrals (optional)
    */ 
    _initialize(atoms_inp, charge);
    s_itg = s_itg_inp;
    t_itg = t_itg_inp;
    v_itg = v_itg_inp;
    e_itg = e_itg_inp;
}    
Molecule::Molecule(int charge, vector<Atom> const& atoms_inp){
    /* Constructor
    Compute integrals implicitly
    Input:
        Variable         Data type                                  Description
        charge_inp       int                                        charge
        atoms_inp        vector<Atom>(n_atom)                       atoms
    */ 
    _initialize(atoms_inp, charge);
    _calcSItg4idx();
    _calcSItg();
    _calcTItg();
    _calcVItg();
    _calcEItg();
}
Molecule::~Molecule(){
    // Destructor
    atoms.clear();
    s_itg.clear();
    t_itg.clear();
    v_itg.clear();
    e_itg.clear();
    s_itg_4idx.clear();
    k_factor.clear();
}
double getSItg(int idx1, int idx2){
    /*
    get overlap integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2)
    double val = s_itg[idx1][idx2];
    return val;
}
double getTItg(int idx1, int idx2){
    /*
    get kinetic energy integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2)
    double val = s_itg[idx1][idx2];
    return val;
}
double getVItg(int idx1, int idx2){
    /*
    get nucleus-electron coulomb integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description
        idx1             int                                        the first index
        idx2             int                                        the second index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2)
    double val = s_itg[idx1][idx2];
    return val;
}    
double getEItg(int idx1, int idx2, int idx3, int idx4){
    /*
    get two-electron integral value between wavefunction idx1 and idx2
    Input:
        Variable         Data type                                  Description            
        idx1             int                                        the first index
        idx2             int                                        the second index
        idx3             int                                        the third index
        idx4             int                                        the forth index
    Return:
        Variable         Data type                                  Description
        val              double                                     integral value
    */
    if(idx2 > idx1) swap(idx1, idx2)
    if(idx4 > idx3) swap(idx3, idx4)
    double val = itg[idx1][idx2][idx3][idx4];
    return val;
}

private:
void _initialize(vector<Atom> const& atoms_inp, int charge){
    // common initialization for explicit and implicit integral
    atoms = atoms_inp;
    n_atom = atoms.size();        
    assert(n_atom > 0 && "Number of atoms should be positive");
    charge = charge_inp;
    _calcNBasis();
    return;
}
void _calcNBasis(){
    // calculate total number of basis functions
    n_basis = 0;
    for(const auto& atom: atoms){
        n_basis += atom.n_basis;
    }
    return;
}
void _calcSItg4idx(){
    // calculate overlap integrals s_itg and 4-indexed overlap integral s_itg_4idx
    is_s_itg_4idx_calculated = true; // set s_itg_4idx as calculated
    for(int atom1_idx = 0 ; atom1_idx < n_atom ; ++atom1){
        
        atom1 = atoms[atom1_idx];
        for(const auto& basis1: atom1.basis_set){
            for(int atom2_idx = 0 ; atom2_idx <= atom1_idx ; ++atom2){
                atom2 = atoms[atom2_idx];
                for(const auto& basis2: atom2.basis_set){

                }
            }
        }
    }
    return;
}
void _checkSItg4idxCalculated(){
    assert(is_s_itg_4idx_calculated && "s_itg_4idx is not calculated, should call _calc_s_itg_4idx() first");
    return;
}
void _calcSItg(){
    // calculate overlap integrals s_itg
    _checkSItg4idxCalculated()
    return;
}
void _calc_TItg(){
    _checkSItg4idxCalculated()
    return;
}
void _calc_VItg(){
    _checkSItg4idxCalculated()
    return;
}
void _calcKFactor(){
    return;
}
void _calcEItg(){
    _calcKFactor();
    return;
}