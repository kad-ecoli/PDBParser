const char* docstring=
"ThreeDblastAlign F1.pdb F2.pdb\n"
"    align two structures by 3D-blast structure alphabet\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NWalign.hpp"
#include "ThreeDblastAlign.hpp"

using namespace std;

int main(int argc, char **argv)
{
    int input_mode=0; // align two fasta files
    int seqID_only=0; // do not just print seqID
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>3)
    {
        input_mode=atoi(argv[3]);
        if (input_mode>=20)
        {
            seqID_only=2;
            input_mode-=20;
        }
    }

    /* parse input */
    vector<string> name_list1,seq_list1;
    vector<string> name_list2,seq_list2;
    vector<vector<int> >seq2int_list1,seq2int_list2; //aa2int
    int seq_num1=read_pdb_as_3dblast(argv[1],name_list1,seq_list1,seq2int_list1);
    int seq_num2=read_pdb_as_3dblast(argv[2],name_list2,seq_list2,seq2int_list2);

    /* do alignment between two files*/
    for (int q=0;q<seq_num1;q++)
    {
        string name1=name_list1[q];
        string seq1=seq_list1[q];
        vector<int> seq2int1=seq2int_list1[q];
        int len1=seq1.length();

        for (int s=0;s<seq_num2;s++)
        {
            string name2=name_list2[s];
            string seq2=seq_list2[s];
            vector<int> seq2int2=seq2int_list2[s];
            int len2=seq2.length();

            string aln1,aln2;
            int aln_score=NWalign(seq1,seq2, seq2int1,seq2int2, aln1,aln2,
                BLOSUM62_3dblast,gapopen_3dblast,gapext_3dblast);

            string aln_str; // colon for identical sequence
            string pos_str; // last digit for position index
            int iden_len,aln_len;  // num of identical/aligned positions
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

            if (seqID_only==2)
            {
                cout<<'>'<<name1<<endl<<aln1<<endl;
                cout<<'>'<<name2<<endl<<aln2<<endl;
                continue;
            }

            cout<<"Length of sequence 1:"<<setw(5)<<len1;
            cout<<" ->"<<name1<<endl;
            cout<<"Length of sequence 2:"<<setw(5)<<len2;
            cout<<" ->"<<name2<<endl;

            cout<<"Alignment score:"<<setw(5)<<aln_score<<endl;
            cout<<"Aligned length:"<<setw(5)<<aln_len<<endl;
            cout<<"Identical length:"<<setw(5)<<iden_len<<endl;
            cout<<"Sequence identity:"<<setiosflags(ios::fixed);
            cout<<setprecision(3)<<setw(9)<<float(iden_len)/len2<<" (=";
            cout<<setw(4)<<iden_len<<'/'<<setw(4)<<len2<<')'<<endl;
            cout<<endl;

            cout<<aln1<<endl;
            cout<<aln_str<<endl;
            cout<<aln2<<endl;
            cout<<pos_str<<endl;

            if (q<seq_num1-1 || s<seq_num2-1) cout<<"$$$$\n"<<endl;
        }
    }
    return 0;
}
