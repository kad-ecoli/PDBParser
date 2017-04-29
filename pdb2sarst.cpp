const char* docstring="pdb2sartst backbone.pdb\n"
"    convert PDB to fasta format SARST code\n";
/* reference:
 * Lo, Wei-Cheng, et al. "Protein structural similarity search by 
 * Ramachandran codes." BMC bioinformatics 8.1 (2007): 307.
 */

#include <iostream>
#include "PDBParser.hpp"
#include "BackboneTorsion.hpp"
#include "FilePathParser.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=1; // only read backbone
    int allowX=1;        // only allow ATOM and MSE

    for (int i=1;i<argc;i++)
    {
        char *infile=argv[i];
        string PDBid=basename_no_ext(infile);
        ModelUnit pdb_entry=read_pdb_structure(infile,atomic_detail,allowX);
        cout<<pdb2sarst(pdb_entry,PDBid);
    }
    return 0;
}
