#include "../class/class.h"

using namespace std;

void setWordMapInput(map<string, string>& );
void setWordMapOutput(map<string, string>& );
Molecule readMolecule(ifstream& , const bool );
double readOutputEnergy(ifstream& );
int readInt(ifstream& );
double readDouble(ifstream& );
void readAtomCoord(ifstream& , vector<Atom>& , const int );
void readBasis(ifstream& , vector<Basis>& , vector<Atom> const& , const int );
void readItg2idx(ifstream& , vector<vector<double>>& , const int);
void readEItg(ifstream& , vector<vector<vector<vector<double>>>>& , const int);