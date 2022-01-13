using namespace std;

void setWordMap(map<string, string>& );
void readMolecule(ifstream& , const bool );
int readInt(ifstream& );
void readAtomCoord(ifstream& , vector<Atom>& , const int );
void readBasis(ifstream& , vector<Basis>& , vector<Atom> const& , const int );
void readItg2idx(ifstream& , vector<vector<double>>& , const int);
void readEItg(ifstream& , vector<vector<vector<vector<double>>>>& , const int);