const char* docstring="reindex_pdb 1 old.pdb new.pdb\n"
"    Renumber residue index in old.pdb such that residue number starts from 1.\n"
"    Output the renumbered pdb to new.pdb\n";

#include <iostream>
#include <cstdlib>
#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    if (argc<3)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=2; // read all atoms
    int allowX=1;        // only allow ATOM and MSE
    int startindex=atoi(argv[1]);    // first residue index in chain

    ModelUnit pdb_entry=read_pdb_structure(argv[2],atomic_detail,allowX);
    reindex_pdb(startindex,pdb_entry);
    write_pdb_structure(argv[3],pdb_entry);
    return 0;
}
