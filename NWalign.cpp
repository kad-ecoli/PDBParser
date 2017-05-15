const char* docstring=""
"Pairwise sequence alignment by standard Needleman-Wunsch algorithm\n"
"    NWalign F1.fasta F2.fasta (align two sequences in fasta file)\n"
"    NWalign F1.pdb F2.pdb 1   (align two sequences in PDB file)\n"
"    NWalign F1.fasta F2.pdb 2 (align Sequence 1 in fasta and 2 in pdb)\n"
"    NWalign GKDGL EVADELVSE 3 (align two sequences typed by keyboard)\n"
"    NWalign GKDGL F.fasta 4   (align sequence 1 by keyboard and 2 in fasta)\n"
"    NWalign GKDGL F.pdb 5     (align sequence 1 by keyboard and 2 in pdb)\n"
"    NWalign F1.fasta ignore 6 (align all sequences within fasta file)\n"
"    NWalign F1.pdb ignore 7   (align all sequences within pdb file)\n"
"\n"
"    NWalign input1 input2 option+10 (only print sequence identity)\n"
"    NWalign input1 input2 option+20 (only print fasta sequence)\n"
"    NWalign input1 input2 option+30 (for each input1 entry, only print\n"
"                                     input2 entry with max identity)\n"
"\n"
//"    NWalign input1 input2 option+100 (glocal-query alignment)\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NWalign.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    int input_mode=0; // align two fasta files
    int glocal=0;     // global or glocal-both alignment
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
        if (input_mode>=100)
        {
            glocal=1;
            input_mode-=100;
        }

        // output format
        if (input_mode>=30)
        {
            seqID_only=3;
            input_mode-=30;
        }
        else if (input_mode>=20)
        {
            seqID_only=2;
            input_mode-=20;
        }
        else if (input_mode>=10)
        {
            seqID_only=1;
            input_mode-=10;
        }
    }

    /* parse input */
    int seq_num1,seq_num2;
    vector<string> name_list1,seq_list1;
    vector<string> name_list2,seq_list2;
    vector<vector<int> >seq2int_list1,seq2int_list2; //aa2int
    switch (input_mode)
    {
        case 0:
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 1:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 2:
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 3:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=get_stdin_seq(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 4:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 5:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2);
            break;
        case 6:
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            break;
        case 7:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,seq_list1,seq2int_list1);
            break;
        default:
            cerr<<"ERROR! Unknown input type "<<input_mode<<endl;
            return 0;
    }

    /* do alignment within one file*/
    if (input_mode==6 || input_mode==7)
    {
        for (int q=0;q<seq_num1-1;q++)
        {
            string name1=name_list1[q];
            string seq1=seq_list1[q];
            vector<int> seq2int1=seq2int_list1[q];
            int len1=seq1.length();

            for (int s=q+1;s<seq_num1;s++)
            {
                string name2=name_list1[s];
                string seq2=seq_list1[s];
                vector<int> seq2int2=seq2int_list1[s];
                int len2=seq2.length();

                string aln1,aln2;
                int aln_score=NWalign(seq1,seq2,seq2int1,seq2int2,aln1,aln2,
                    BLOSUM62,gapopen_blosum62,gapext_blosum62,glocal);

                string aln_str; // colon for identical sequence
                string pos_str; // last digit for position index
                int iden_len,aln_len;  // num of identical/aligned positions
                get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

                if (seqID_only==2) // fasta alignment only
                {
                    cout<<'>'<<name1<<endl<<aln1<<endl;
                    cout<<'>'<<name2<<endl<<aln2<<endl;
                    if (q<seq_num1-1 || s<seq_num2-1) cout<<"$$$$\n"<<endl;
                }
                else if (seqID_only==1) // seqID only
                {
                    cout<<name1<<'\t'<<name2<<'\t';
                    cout<<setiosflags(ios::fixed)<<setprecision(4);
                    cout<<float(iden_len)/len1<<'\t';
                    cout<<setiosflags(ios::fixed)<<setprecision(4);
                    cout<<float(iden_len)/len2<<endl;
                }
                else // full output
                {
                    cout<<"Length of sequence 1:"<<setw(5)<<len1;
                    cout<<" ->"<<name1<<endl;
                    cout<<"Length of sequence 2:"<<setw(5)<<len2;
                    cout<<" ->"<<name2<<endl;

                    cout<<"Alignment score:"<<aln_score<<endl;
                    cout<<"Aligned length:"<<setw(5)<<aln_len<<endl;
                    cout<<"Identical length:"<<setw(5)<<iden_len<<endl;
                    cout<<"Sequence identity:"<<setiosflags(ios::fixed);
                    cout<<setprecision(3)<<setw(9)<<float(iden_len)/len2;
                    cout<<" (="<<setw(4)<<iden_len<<'/'<<setw(4)<<len2;
                    cout<<")\n\n";

                    cout<<aln1<<endl;
                    cout<<aln_str<<endl;
                    cout<<aln2<<endl;
                    cout<<pos_str<<endl;

                    if (q<seq_num1-2 || s<seq_num1-1) cout<<"$$$$\n"<<endl;
                }
            }
        }
        return 0;
    }

    /* do alignment between two files*/
    for (int q=0;q<seq_num1;q++)
    {
        string name1=name_list1[q];
        string seq1=seq_list1[q];
        vector<int> seq2int1=seq2int_list1[q];
        int len1=seq1.length();
        int max_iden_len=0; // identical positions with max seqID seq2
        string max_seqID_name2="";
        int max_seqID_len2=0;

        for (int s=0;s<seq_num2;s++)
        {
            string name2=name_list2[s];
            string seq2=seq_list2[s];
            vector<int> seq2int2=seq2int_list2[s];
            int len2=seq2.length();

            string aln1,aln2;
            int aln_score=NWalign(seq1,seq2, seq2int1,seq2int2, aln1,aln2,
                BLOSUM62,gapopen_blosum62,gapext_blosum62,glocal);

            string aln_str; // colon for identical sequence
            string pos_str; // last digit for position index
            int iden_len,aln_len;  // num of identical/aligned positions
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

            if (seqID_only==3) // max seqID only
            {
                if (max_iden_len<=iden_len)
                {
                    max_iden_len=iden_len;
                    max_seqID_name2=name2;
                    max_seqID_len2=len2;
                }
            }
            else if (seqID_only==2) // fasta alignment only
            {
                cout<<'>'<<name1<<endl<<aln1<<endl;
                cout<<'>'<<name2<<endl<<aln2<<endl;
                if (q<seq_num1-1 || s<seq_num2-1) cout<<"$$$$\n"<<endl;
            }
            else if (seqID_only==1) // seqID only
            {
                cout<<name1<<'\t'<<name2<<'\t';
                cout<<setiosflags(ios::fixed)<<setprecision(4);
                cout<<float(iden_len)/len1<<'\t';
                cout<<setiosflags(ios::fixed)<<setprecision(4);
                cout<<float(iden_len)/len2<<endl;
                if (seq_num1>1 && s==seq_num2-1) cout<<"$$$$\n";
            }
            else // full output
            {
                cout<<"Length of sequence 1:"<<setw(5)<<len1;
                cout<<" ->"<<name1<<endl;
                cout<<"Length of sequence 2:"<<setw(5)<<len2;
                cout<<" ->"<<name2<<endl;

                cout<<"Alignment score:"<<setw(5)<<aln_score<<endl;
                cout<<"Aligned length:"<<setw(5)<<aln_len<<endl;
                cout<<"Identical length:"<<setw(5)<<iden_len<<endl;
                cout<<"Sequence identity:"<<setiosflags(ios::fixed);
                cout<<setprecision(3)<<setw(9)<<float(iden_len)/len2;
                cout<<" (="<<setw(4)<<iden_len<<'/'<<setw(4)<<len2;                
                cout<<")\n\n";

                cout<<aln1<<endl;
                cout<<aln_str<<endl;
                cout<<aln2<<endl;
                cout<<pos_str<<endl;
                if (q<seq_num1-1 || s<seq_num2-1) cout<<"$$$$\n"<<endl;
            }

        }
        
        if (seqID_only==3)
        {
            cout<<name1<<'\t'<<max_seqID_name2<<'\t';
            cout<<setiosflags(ios::fixed)<<setprecision(4);
            cout<<float(max_iden_len)/len1<<'\t';
            cout<<setiosflags(ios::fixed)<<setprecision(4);
            cout<<float(max_iden_len)/max_seqID_len2<<endl;
        }
    }
    return 0;
}
