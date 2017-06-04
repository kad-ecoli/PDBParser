const char* docstring=""
"Pairwise sequence alignment by standard Needleman-Wunsch algorithm\n"
"    NWalign F1.fasta F2.fasta (align two sequences in fasta file)\n"
"    NWalign F1.pdb F2.pdb   1 (align two sequences in PDB file)\n"
"    NWalign F1.fasta F2.pdb 2 (align Sequence 1 in fasta and 2 in pdb)\n"
"    NWalign GKDGL EVADELVSE 3 (align two sequences typed by keyboard)\n"
"    NWalign GKDGL F.fasta   4 (align sequence 1 by keyboard and 2 in fasta)\n"
"    NWalign GKDGL F.pdb     5 (align sequence 1 by keyboard and 2 in pdb)\n"
"    NWalign F1.fasta ignore 6 (align all sequences within fasta file)\n"
"    NWalign F1.pdb ignore   7 (align all sequences within pdb file)\n"
"\n"
"    NWalign input1 input2 option+10 (only print sequence identity)\n"
"    NWalign input1 input2 option+20 (only print fasta sequence)\n"
"    NWalign input1 input2 option+30 (for each input1 entry, only print\n"
"                                     input2 entry with max identity)\n"
"\n"
"    NWalign input1 input2 option+100 (glocal-query alignment)\n"
"    NWalign input1 input2 option+200 (glocal-both alignment)\n"
"    NWalign input1 input2 option+300 (local alignment)\n"
"\n"
"    NWalign input1 input2 option+1000 (using chi-1 rotamer sequence)\n"
"    NWalign input1 input2 option+2000 (using SARST sequence)\n"
"    NWalign input1 input2 option+3000 (using 3d-blast sequence)\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NWalign.hpp"
#include "ROTSUMalign.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    int input_mode=0; // align two fasta files
    int glocal=0;     // global or glocal-both alignment
    int seqID_only=0; // do not just print seqID
    int seq_type=0;   // amino acid sequence
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>3)
    {
        input_mode=atoi(argv[3]);
        // sequence type
        if (input_mode>=1000)
        {
            seq_type=int(input_mode/1000);
            input_mode=input_mode % 1000;
        }

        // alignment algorithm
        if (input_mode>=100)
        {
            glocal=int(input_mode/100);
            input_mode=input_mode % 100;
        }

        // output format
        if (input_mode>=10)
        {
            seqID_only=int(input_mode/10);
            input_mode=input_mode % 10;
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
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=read_fasta(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 1:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 2:
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 3:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=get_stdin_seq(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 4:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=read_fasta(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 5:
            seq_num1=get_stdin_seq(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,seq_list2,seq2int_list2,seq_type);
            break;
        case 6:
            seq_num1=read_fasta(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            break;
        case 7:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,seq_list1,seq2int_list1,seq_type);
            break;
        default:
            cerr<<"ERROR! Unknown input type "<<input_mode<<endl;
            return 0;
    }
    if (input_mode>=6) seq_num2=seq_num1;

    /* matrix for all-against-all seqID within one file */
    vector<vector<int> > iden_len_all_against_all;
    vector<vector<int> > aln_len_all_against_all;
    if (input_mode>=6 && seqID_only==3 && glocal!=1)
    {
        vector<int> temp_int(seq_num2,1);
        iden_len_all_against_all.assign(seq_num1,temp_int);
        aln_len_all_against_all.assign(seq_num1,temp_int);
        temp_int.clear();
    }

    /* do alignment */
    int q,s,len1,len2,max_iden_len,max_aln_len,max_seqID_len2;
    string name1,name2,seq1,seq2,max_seqID_name2;
    vector<int> seq2int1,seq2int2;
    for (q=0;q<((input_mode>=6 && glocal!=1)?(seq_num1-1):seq_num1);q++)
    {
        name1=name_list1[q];
        seq1=seq_list1[q];
        seq2int1=seq2int_list1[q];
        len1=seq1.length();
        max_iden_len=0; // identical positions with max seqID seq2
        max_aln_len=0;  // aligned   positions with max seqID seq2
        max_seqID_name2="";
        max_seqID_len2=0;

        if (q>0 && input_mode>=6 && seqID_only==3 && glocal!=1)
        {   // find max seqID pair from old alignment
            for (s=0;s<q;s++)
            {
                if (iden_len_all_against_all[s][q]>max_aln_len)
                {
                    max_iden_len=iden_len_all_against_all[s][q];
                    max_aln_len=aln_len_all_against_all[s][q];
                    max_seqID_name2=name_list1[s];
                    max_seqID_len2=seq_list1[s].length();
                }
            }
        }

        for (s=((input_mode>=6 && glocal!=1)?(q+1):0);s<seq_num2;s++)
        {
            if (s==q && input_mode>=6) continue; // glocal==1
            if (input_mode>=6)
            {
                name2=name_list1[s];
                seq2=seq_list1[s];
                seq2int2=seq2int_list1[s];
            }
            else
            {
                name2=name_list2[s];
                seq2=seq_list2[s];
                seq2int2=seq2int_list2[s];
            }
            len2=seq2.length();

            string aln1,aln2;
            int aln_score;
            switch (seq_type)
            {
                case 1: // chi-1 rotamer
                    aln_score=NWalign(seq1,seq2, seq2int1,seq2int2,aln1,
                        aln2, ROTSUM8,gapopen_rotsum8,gapext_rotsum8,glocal);
                    break;
                case 2: // sarst code
                    aln_score=NWalign(seq1,seq2,seq2int1,seq2int2,aln1,aln2,
                        BLOSUM62_sarst,gapopen_sarst,gapext_sarst,glocal);
                    break;
                case 3: // 3d-blast sequence
                    aln_score=NWalign(seq1,seq2,seq2int1,seq2int2,aln1,aln2,
                        BLOSUM62_3dblast,gapopen_3dblast,gapext_3dblast,
                        glocal);
                    break;
                default: // amino acid
                    aln_score=NWalign(seq1,seq2,seq2int1,seq2int2,aln1,aln2,
                        BLOSUM62,gapopen_blosum62,gapext_blosum62,glocal);
            }

            string aln_str; // colon for identical sequence
            string pos_str; // last digit for position index
            int iden_len,aln_len;  // num of identical/aligned positions
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

            if (seqID_only==3) // max seqID only
            {
                if (input_mode>=6 && glocal!=1)
                {
                    iden_len_all_against_all[q][s]=
                    iden_len_all_against_all[s][q]=iden_len;
                    aln_len_all_against_all[q][s]=
                    aln_len_all_against_all[s][q]=aln_len;
                }
                if (max_iden_len<=iden_len)
                {
                    max_iden_len=iden_len;
                    max_aln_len=aln_len;
                    max_seqID_name2=name2;
                    max_seqID_len2=len2;
                }
            }
            else if (seqID_only==2) // fasta alignment only
            {
                cout<<'>'<<name1<<endl<<aln1<<endl;
                cout<<'>'<<name2<<endl<<aln2<<endl;
                if (((input_mode>=6)?(q<seq_num1-2):(q<seq_num1-1)) || 
                    s<seq_num2-1) cout<<"$$$$\n"<<endl;
            }
            else if (seqID_only==1) // seqID only
            {
                cout<<name1<<'\t'<<name2<<'\t'
                    <<setiosflags(ios::fixed)<<setprecision(4)
                    <<float(iden_len)/len1<<'\t'
                    <<float(iden_len)/len2<<'\t'
                    <<float(iden_len)/aln_len<<endl;
                if (input_mode<6 && q<seq_num1-1 && s==seq_num2-1)
                    cout<<"$$$$\n";
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
                if (((input_mode>=6)?(q<seq_num1-2):(q<seq_num1-1)) || 
                    s<seq_num2-1) cout<<"$$$$\n"<<endl;
            }

        }
        
        if (seqID_only==3)
        {
            cout<<name1<<'\t'<<max_seqID_name2<<'\t'
                <<setiosflags(ios::fixed)<<setprecision(4)
                <<float(max_iden_len)/len1<<'\t'
                <<float(max_iden_len)/max_seqID_len2<<'\t'
                <<float(max_iden_len)/max_aln_len<<endl;
        }
    }
    return 0;
}
