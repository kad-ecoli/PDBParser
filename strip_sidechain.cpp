const char* docstring="strip_sidechain full_atom.pdb backbone.pdb\n"
"    remove sidechain and hydrogen atoms from full_atom.pdb,\n"
"    saving CA, C, N, O atoms to backbone.pdb\n";

#include <iostream>
#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    if (argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=1; // only read backbone
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);
    write_pdb_structure(argv[2],pdb_entry);
    return 0;
}
