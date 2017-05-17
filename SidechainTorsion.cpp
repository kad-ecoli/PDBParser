const char* docstring="SidechainTorsion full.pdb\n"
"    calculate sidechain torsion angles (chi1, chi2, chi3, chi4) from full.pdb\n";

#include <iostream>
#include "PDBParser.hpp"
#include "SidechainTorsion.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    
    int atomic_detail=2; // only read full atom structure
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argv[1],atomic_detail,allowX);

    int c,r; // chain index, residue index;
    cout<<" AA c resi     chi1    chi2    chi3    chi4"<<endl;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        vector<vector<float> >chi_mat=SidechainTorsion(pdb_entry.chains[c]);
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            cout<<setw(3)<<pdb_entry.chains[c].residues[r].resn<<' '
                <<pdb_entry.chains[c].chainID<<' '
                <<setw(4)<<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode;
            if (chi_mat[r][0]==360)
            {
                cout<<endl;
                continue;
            }
            cout<<' '<<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<chi_mat[r][0];
            if (chi_mat[r][1]==360)
            {
                cout<<endl;
                continue;
            }
            cout<<' '<<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<chi_mat[r][1];
            if (chi_mat[r][2]==360)
            {
                cout<<endl;
                continue;
            }
            cout<<' '<<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<chi_mat[r][2];
            if (chi_mat[r][3]==360)
            {
                cout<<endl;
                continue;
            }
            cout<<' '<<setiosflags(ios::fixed)<<setprecision(2)
                <<setw(7)<<chi_mat[r][3]<<endl;
        }
        chi_mat.clear();
    }
    return 0;
}
