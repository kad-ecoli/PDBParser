const char* docstring="backbone_torsion backbone.pdb\n"
"    calculate backbone torsion angles (Omega, Phi, Psi) from backbone.pdb\n";

#include <iostream>
#include "PDBParser.hpp"
#include "MathTools.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=1; // only read backbone
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);

    int c,r; // chain index, residue index;
    cout<<" AA c resi      omg     phi     psi"<<endl;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        vector<vector<float> > angle_mat=BackBoneTorsion(pdb_entry.chains[c]);
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            cout<<setw(3)<<pdb_entry.chains[c].residues[r].resn<<' '
                <<pdb_entry.chains[c].chainID<<' '
                <<setw(4)<<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode<<' '
                <<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<angle_mat[r][0]<<' '
                <<setw(7)<<angle_mat[r][1]<<' '
                <<setw(7)<<angle_mat[r][2]<<endl;
        }
    }
    return 0;
}
