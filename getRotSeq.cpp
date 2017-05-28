const char* docstring="getRotSeq full.pdb\n"
"    get fasta format rotamer sequence based on chi1 angle\n";

#include <iostream>
#include "getRotSeq.hpp"
#include "FilePathParser.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=2; // read full atom structure
    int allowX=1;        // only allow ATOM and MSE

    for (int i=1;i<argc;i++)
    {
        ModelUnit pdb_entry=read_pdb_structure(argv[i],atomic_detail,allowX);
        string PDBid=basename_no_ext(argv[i]);
        for (int c=0;c<pdb_entry.chains.size();c++)
        {
            cout<<'>'+PDBid+':'+pdb_entry.chains[c].chainID_full<<endl;
            cout<<getRotSeq(pdb_entry.chains[c])<<endl;
        }
    }
    return 0;
}
