const char* docstring="pdb2fasta pdb.pdb > seq.fasta\n"
"    convert PDB file 'pdb.pdd' to fasta file 'seq.fasta'\n";

#include <iostream>
#include "PDBParser.hpp"
#include "FilePathParser.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=0; // only read CA
    int allowX=1;        // only allow ATOM and MSE

    for (int i=1;i<argc;i++)
    {
        char *infile=argv[i];
        string PDBid=basename_no_ext(infile);
        ModelUnit pdb_entry=read_pdb_structure(infile,atomic_detail,allowX);
        string txt=pdb2fasta(pdb_entry,PDBid);
        cout<<txt;
    }
    return 0;
}
