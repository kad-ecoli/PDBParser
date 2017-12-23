const char* docstring=""
"ccealign: Combinatorial Extension (CE) method for pairwise structure alignment\n"
"    ccealign F1.pdb F2.pdb\n"
"\n"
"options:\n"
"    -prefix=CE_ - output superposed PDB, using 'CE_' as file prefix\n"
;

#include "PDBParser.hpp"
#include "FilePathParser.hpp"

#include "ccealign.hpp"
#include "ce2tm.hpp"

using namespace std;

int main(int argc, char **argv)
{
    string prefix=""; // do not output superposed PDB by default
    vector<string> argv_list; // list of PDB input
    int arg,i;
    for (arg=1;arg<argc;arg++)
    {
        if (string(argv[arg]).substr(0,8)=="-prefix=")
            prefix=string(argv[arg]).substr(8);
        else
            argv_list.push_back(argv[arg]);
    }
    if (argv_list.size()<=1)
    {
        cout<<docstring<<endl;
        return 0;
    }

    vector<ChainUnit> chain_list;
    for (i=0;i<argv_list.size();i++)
    {
        chain_list.push_back(  // 0 - CA only, 1 - ATOM & MSE
            read_pdb_structure(argv_list[i].c_str(),0,1).chains[0]);
        pdb2fasta(chain_list[i]);
    }

    string aln_moving, aln_fixed, aln_str;
    float rmsd, tmscore_moving, tmscore_fixed, seqID;
    int L_ali;
    for (i=0;i<chain_list.size()-1;i++)
    {
        if (i>=1) cout<<"$$$$\n\n";

        vector<vector<float> > rVal=ccealign_ccealign(
            chain_list[i], chain_list[chain_list.size()-1]);

        cout<<"Name of Chain_1: "<<argv_list[i]<<endl
            <<"Name of Chain_2: "<<argv_list[argv_list.size()-1]<<endl
            <<"Length of Chain_1: "
            <<chain_list[i].residues.size()<<" residues"<<endl
            <<"Length of Chain_2: "
            <<chain_list[chain_list.size()-1].residues.size()<<" residues"<<endl
            <<endl;

        super2aln(chain_list[i], chain_list[chain_list.size()-1],
            rVal, aln_moving, aln_fixed, aln_str, 
            rmsd, tmscore_moving, tmscore_fixed, L_ali, seqID);

        cout<<"Aligned length= "<<L_ali<<", RMSD= "
            <<setiosflags(ios::fixed)<<setprecision(2)<<rmsd
            <<", Seq_ID=n_identical/n_aligned= "<<setprecision(3)<<seqID
            <<endl<<setprecision(5)
            <<"TM-score= "<<tmscore_moving
            <<" (if normalized by length of Chain_1)"<<endl
            <<"TM-score= "<<tmscore_fixed
            <<" (if normalized by length of Chain_2)"<<endl
            <<endl
            <<"(':' denotes aligned residue pairs of d < 5.0 A,"
            <<" '.' denotes other aligned residues)"<<endl
            <<aln_moving<<endl<<aln_str<<endl<<aln_fixed<<endl;

        if (prefix!="")
        {
            string outfile=prefix+basename(argv_list[i].c_str());
            ChainUnit old_chain=read_pdb_structure(argv_list[i].c_str(),2,3
                ).chains[0]; // 2 - all atoms; 3 - no residue conversion
            ChainUnit super_chain=old_chain;
            int r,a;

            for (r=0;r<super_chain.residues.size();r++)
            {
                for (a=0;a<super_chain.residues[r].atoms.size();a++)
                {
                    transform_selection(old_chain.residues[r].atoms[a].xyz,
                        rVal, super_chain.residues[r].atoms[a].xyz);
                }
            }
            write_pdb_structure(outfile.c_str(),super_chain);

            super_chain.residues.clear();
        }
    }
    return 0;
}
