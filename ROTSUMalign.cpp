const char* docstring=
"pairwise structure alignment by their chi-1 rotamer sequence\n"
"    ROTSUMalign F1.fasta F2.fasta (align two rotamer sequences in fasta file)\n"
"    ROTSUMalign F1.pdb   F2.pdb 1 (align two structures in PDB file)\n"
"    ROTSUMalign F1.fasta F2.pdb 2 (align rotamer sequence 1 in fasta and 2 in pdb)\n"
"    ROTSUMalign GKDGL EVADELVSE 3 (align two rotamer sequences typed by keyboard)\n"
"    ROTSUMalign GKDGL   F.fasta 4 (align rotamer sequence 1 by keyboard and 2 in fasta)\n"
"    ROTSUMalign GKDGL   F.pdb   5 (align rotamer sequence 1 by keyboard and 2 in pdb)\n"
"\n"
"    ROTSUMalign input1 input2 option+20 (only print fasta sequence)\n"
"\n"
"    ROTSUMalign input1 input2 option+100 (glocal-query alignment)\n"
"    ROTSUMalign input1 input2 option+200 (glocal-both alignment)\n"
"    ROTSUMalign input1 input2 option+300 (local alignment)\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "ROTSUMalign.hpp"

using namespace std;

int main(int argc, char **argv)
{
    int input_mode=0; // align two fasta files
    int glocal=0;     // global, glocal, or local alignment
    int seqID_only=0; // do not just print seqID
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>3)
    {
        input_mode=atoi(argv[3]);
        // alignment algorithm
        if (input_mode>=300)
        {
            glocal=3;
            input_mode-=300;
        }
        else if (input_mode>=200)
        {
            glocal=2;
            input_mode-=200;
        }
        else if (input_mode>=100)
        {
            glocal=1;
            input_mode-=100;
        }

        // output format
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
    int seq_num1=0;
    int seq_num2=0;

    switch (input_mode)
    {
        case 0:
            seq_num1=read_rotseq_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_rotseq_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 1:
            seq_num1=read_pdb_as_rotseq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_rotseq(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 2:
            seq_num1=read_rotseq_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_rotseq(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 3:
            seq_num1=get_stdin_rotseq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=get_stdin_rotseq(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 4:
            seq_num1=get_stdin_rotseq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_rotseq_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 5:
            seq_num1=get_stdin_rotseq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_rotseq(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        default:
            cerr<<"ERROR! Unknown input type "<<input_mode<<endl;
            return 0;
    }

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
            float aln_score=NWalign(seq1,seq2, seq2int1,seq2int2, aln1,aln2,
                ROTSUM8,gapopen_rotsum8,gapext_rotsum8,glocal);

            string aln_str; // colon for identical sequence
            string pos_str; // last digit for position index
            int iden_len,aln_len;  // num of identical/aligned positions
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

            if (seqID_only==2)
            {
                cout<<'>'<<name1<<endl<<aln1<<endl;
                cout<<'>'<<name2<<endl<<aln2<<endl;
                if (q<seq_num1-1 || s<seq_num2-1) cout<<"$$$$\n"<<endl;
            }
            else
            {
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
    }
    return 0;
}
