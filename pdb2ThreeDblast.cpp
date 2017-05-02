const char* docstring="pdb2ThreeDblast CA.pdb\n"
"    convert PDB to fasta format structural alphabet in 3D-blast\n";
/* reference:
 * Tung, Chi-Hua, Jhang-Wei Huang, and Jinn-Moon Yang. "Kappa-alpha plot
 * derived structural alphabet and BLOSUM-like substitution matrix for rapid
 * search of protein structure database." Genome biology 8.3 (2007): R31.
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
        cout<<pdb2ThreeDblast(pdb_entry,PDBid);
    }
    return 0;
}
