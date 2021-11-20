const char* docstring="centerModel old.pdb new.pdb\n"
"    translate the structure so that the center is (0,0,0)\n"
;

#include <iostream>
#include <cstdlib>
#include "PDBParser.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=2; // read all atoms
    int allowX=3;        // all residues, no conversion
    
    int c,r,a;
    int atom_count=0;
    double x0=0;
    double y0=0;
    double z0=0;
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
            {
                x0+=pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
                y0+=pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
                z0+=pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
                atom_count++;
            }
        }
    }
    x0/=atom_count;
    y0/=atom_count;
    z0/=atom_count;
    x0=int(x0*1000+0.5)/1000.; // only keep 3 digits after decimal
    y0=int(y0*1000+0.5)/1000.; // only keep 3 digits after decimal
    z0=int(z0*1000+0.5)/1000.; // only keep 3 digits after decimal
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
            {
                pdb_entry.chains[c].residues[r].atoms[a].xyz[0]-=x0;
                pdb_entry.chains[c].residues[r].atoms[a].xyz[1]-=y0;
                pdb_entry.chains[c].residues[r].atoms[a].xyz[2]-=z0;
            }
        }
    }
    write_pdb_structure((argc>2?argv[2]:"-"),pdb_entry);
    vector<ChainUnit>().swap(pdb_entry.chains);
    return 0;
}
