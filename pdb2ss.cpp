const char* docstring="pdb2ss CA.pdb\n"
"    convert PDB to fasta format secondary structure\n";
/* reference:
 * Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)
 */

#include <iostream>
#include "PDBParser.hpp"
#include "StructuralAlphabet.hpp"
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
        cout<<pdb2ss(pdb_entry,PDBid);
    }
    return 0;
}
