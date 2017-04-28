const char* docstring="split_chain pdb.pdb\n"
"    split multichain PDB into individual chains.\n";

#include <iostream>
#include <cstdlib>
#include "PDBParser.hpp"
#include "FilePathParser.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=2; // read all atoms
    int allowX=3;        // all residues, no conversion

    for (int i=1;i<argc;i++)
    {
        ModelUnit pdb_entry=read_pdb_structure(argv[i],atomic_detail,allowX);

        string filename_prefix=filename_no_ext(argv[i]);
        for (int c=0;c<pdb_entry.chains.size();c++)
        {
            string filename=filename_prefix+
                pdb_entry.chains[c].chainID_full+".pdb";
            cout<<filename<<endl;
            write_pdb_structure(filename.c_str(), pdb_entry.chains[c]);
        }
    }
    return 0;
}
