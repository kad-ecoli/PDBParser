const char* docstring=
"RMSD superposition of two proteins aligned by sequence\n"
"    pdb2rmsd F1.pdb F2.pdb\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NWalign.hpp"
#include "PDBParser.hpp"
#include "pdb2rmsd.hpp"
#include "Superpose.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=0; // only read CA
    int allowX=1;        // only allow ATOM and MSE

    /* parse input */
    string PDBid1=filename_no_ext(argv[1]);
    string PDBid2=filename_no_ext(argv[2]);
    ModelUnit pdb_entry1=read_pdb_structure(argv[1],atomic_detail,allowX);
    ModelUnit pdb_entry2=read_pdb_structure(argv[2],atomic_detail,allowX);

    /* pdb2fasta */
    int seq_num1=pdb_entry1.chains.size();
    int seq_num2=pdb_entry2.chains.size();
    vector<vector<int> > seq2int1_list;
    vector<vector<int> > seq2int2_list;
    for (int c=0;c<seq_num1;c++)
    {
        pdb2fasta(pdb_entry1.chains[c]);
        seq2int1_list.push_back(aa2int(pdb_entry1.chains[c].sequence));
    }
    for (int c=0;c<seq_num2;c++)
    {
        pdb2fasta(pdb_entry2.chains[c]);
        seq2int2_list.push_back(aa2int(pdb_entry2.chains[c].sequence));
    }

    /* NWalign + RMSD*/
    int c1,c2; // index of chain
    string aln1,aln2; // alignment
    int L1,L2; // sequence length
    int aln_len; // aligned position number
    vector<double> tmp_array(3,0.);
    vector<vector<double> > xyz_list1,xyz_list2; // coordinate of aligned residue
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    int r; // resi
    string aln_str; // string to indicate match/mistmatch
    string pos_str; // string to indicate position
    int iden_len; // number of identical position

    for (c1=0;c1<seq_num1;c1++)
    {
        for (c2=0;c2<seq_num2;c2++)
        {
            /* NWalign */
            NWalign(pdb_entry1.chains[c1].sequence,
                pdb_entry2.chains[c2].sequence,
                seq2int1_list[c1],seq2int2_list[c2], aln1,aln2);
            cout<<'>'<<PDBid1<<':'<<pdb_entry1.chains[c1].chainID_full
                <<'\t'<<pdb_entry1.chains[c1].sequence.length()<<endl
                <<aln1<<endl;
            cout<<'>'<<PDBid2<<':'<<pdb_entry2.chains[c2].chainID_full
                <<'\t'<<pdb_entry2.chains[c2].sequence.length()<<endl
                <<aln2<<endl;
            
            /* seqID */
            L1=pdb_entry1.chains[c1].sequence.length();
            L2=pdb_entry2.chains[c2].sequence.length();
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);
            aln_str.clear();
            pos_str.clear();
            cout<<setiosflags(ios::fixed)<<setprecision(4)
                <<"# IDali="<<1.*iden_len/aln_len
                <<"\tidentity1="<<1.*iden_len/L1
                <<"\tidentity2="<<1.*iden_len/L2<<endl;
            
            /* extract xyz coordinate */
            aln2coor(aln1,aln2,pdb_entry1.chains[c1],pdb_entry2.chains[c2],
                xyz_list1,xyz_list2,atomic_detail);
            cout<<setiosflags(ios::fixed)<<setprecision(4)
                <<"# Lali="<<aln_len<<"\tcoverage1="<<1.*aln_len/L1
                <<"\tcoverage2="<<1.*aln_len/L2<<endl;

            /* RMSD superposition */
            if (aln_len!=0)
            {
                /* kabsch */
                RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

                /* change coordinate */
                vector<vector<double> > super_xyz_list1(aln_len,tmp_array);
                for(r=0; r<aln_len; r++)
	                ChangeCoor(xyz_list1[r], RotMatix, TranVect, 
                        super_xyz_list1[r]);

                /* RMSD */
                cout<<setiosflags(ios::fixed)<<setprecision(4)
                    <<"# RMSD="<<calRMSD(super_xyz_list1, xyz_list2)
                    <<"\tTM-score1="<<calTMscore(super_xyz_list1,xyz_list2,L1)
                    <<"\tTM-score2="<<calTMscore(super_xyz_list1,xyz_list2,L2)
                    <<endl;
                
                /* clean up */
                RotMatix.clear();
                TranVect.clear();
                super_xyz_list1.clear();
            }

            /* clean up */
            xyz_list1.clear();
            xyz_list2.clear();
            aln1.clear();
            aln2.clear();

            if (c1<seq_num1-1 || c2<seq_num2-1) cout<<"$$$$\n"<<endl;
        }
    }

    /* clean up*/
    seq2int1_list.clear();
    seq2int2_list.clear();
    pdb_entry1.chains.clear();
    pdb_entry2.chains.clear();
    PDBid1.clear();
    PDBid2.clear();
    return 0;
}
