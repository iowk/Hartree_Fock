#include <bits/stdc++.h>
#include "read_input.h"
#include "../class/class.h"

using namespace std;

int main(){
    string inp_path;
    cout << "Enter input file path:" << endl;
    getline(cin, inp_path);
    ifstream inp_file(inp_path);
    if(inp_file.is_open()){
        string s_read_itg;
        bool is_read_itg = false;
        cout << "Read integrals from input file (Y/n)?" << endl;
        getline(cin, s_read_itg);
        if(s_read_itg.size() == 0 || s_read_itg[0] == 'Y' || s_read_itg[0] == 'y'){
            is_read_itg = true;
            cout << "Integrals will be read from input file explicitly" << endl;
        }
        else{
            cout << "Integrals will be calculated implicitly" << endl;
        }
        readMolecule(inp_file, is_read_itg);
    }
    else{
        cout << "Unable to open file" << endl;
        return 0;
    }    
}