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
"    NWalign input1 input2 option+4000 (using secondary structure by TMalign)\n"
"\n"
"    NWalign input1 input2 option+10000 (RMSD superposition)\n"
"    NWalign input1 input2 option+20000 (TM-score superposition)\n"
"    NWalign input1 input2 option+30000 (fast TM-score superposition)\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "NWalign.hpp"
#include "ROTSUMalign.hpp"
#include "pdb2rmsd.hpp"
#include "Superpose.hpp"
#include "TMalign.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    int input_mode=0; // align two fasta files
    int glocal=0;     // global or glocal-both alignment
    int seqID_only=0; // do not just print seqID
    int seq_type=0;   // amino acid sequence
    int super_type=0; // do not perform superposition
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>3)
    {
        input_mode=atoi(argv[3]);
        // superposition type
        if (input_mode>=10000)
        {
            super_type=int(input_mode/10000);
            input_mode=    input_mode%10000;
        }

        // sequence type
        if (input_mode>=1000)
        {
            seq_type=int(input_mode/1000);
            input_mode=  input_mode%1000;
        }

        // alignment algorithm
        if (input_mode>=100)
        {
            glocal=int(input_mode/100);
            input_mode=input_mode%100;
        }

        // output format
        if (input_mode>=10)
        {
            seqID_only=int(input_mode/10);
            input_mode=    input_mode%10;
        }
    }

    /* check if input is superposable */
    if (super_type>0 && input_mode!=1 && input_mode!=7)
    {
        cerr<<"ERROR! Cannot superpose input format type "<<input_mode<<endl;
        return 0;
    }

    /* parse input */
    int seq_num1,seq_num2;
    vector<string> name_list1,seq_list1;
    vector<string> name_list2,seq_list2;
    vector<vector<int> >seq2int_list1,seq2int_list2; //aa2int
    ModelUnit pdb_entry1,pdb_entry2;
    switch (input_mode)
    {
        case 0:
            seq_num1=read_fasta(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            seq_num2=read_fasta(argv[2],name_list2,
                seq_list2,seq2int_list2,seq_type);
            break;
        case 1:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,
                seq_list1,seq2int_list1,pdb_entry1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,
                seq_list2,seq2int_list2,pdb_entry2,seq_type);
            break;
        case 2:
            seq_num1=read_fasta(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,
                seq_list2,seq2int_list2,pdb_entry2,seq_type);
            break;
        case 3:
            seq_num1=get_stdin_seq(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            seq_num2=get_stdin_seq(argv[2],name_list2,
                seq_list2,seq2int_list2,seq_type);
            break;
        case 4:
            seq_num1=get_stdin_seq(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            seq_num2=read_fasta(argv[2],name_list2,
                seq_list2,seq2int_list2,seq_type);
            break;
        case 5:
            seq_num1=get_stdin_seq(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            seq_num2=read_pdb_as_fasta(argv[2],name_list2,
                seq_list2,seq2int_list2,pdb_entry2,seq_type);
            break;
        case 6:
            seq_num1=read_fasta(argv[1],name_list1,
                seq_list1,seq2int_list1,seq_type);
            break;
        case 7:
            seq_num1=read_pdb_as_fasta(argv[1],name_list1,
                seq_list1,seq2int_list1,pdb_entry1,seq_type);
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

    /* variables for RMSD */
    vector<float> tmp_array(3,0.);
    vector<vector<float> > xyz_list1,xyz_list2; // coordinate of aligned residue
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t
    float rmsd,tmscore1,tmscore2;

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
            int aln_score=NWalign(seq1,seq2, seq2int1,seq2int2, aln1,aln2,
                seq_type,glocal);

            string aln_str; // colon for identical sequence
            string pos_str; // last digit for position index
            int iden_len,aln_len;  // num of identical/aligned positions
            get_seqID(aln1,aln2,aln_str,pos_str,iden_len,aln_len);

            /* RMSD/TM-score superposition */
            if (aln_len==0 && super_type>=1)
            {
                rmsd=0;
                tmscore1=0;
                tmscore2=0;
            }
            else if (aln_len!=0 && super_type==1) // RMSD
            {
                if (input_mode>=6)
                    aln2coor(aln1,aln2,pdb_entry1.chains[q],pdb_entry1.chains[s],
                        xyz_list1,xyz_list2,1);
                else
                    aln2coor(aln1,aln2,pdb_entry1.chains[q],pdb_entry2.chains[s],
                        xyz_list1,xyz_list2,1);

                /* kabsch */
                RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

                /* change coordinate */
                vector<vector<float> > super_xyz_list1(aln_len,tmp_array);
                for(int r=0; r<aln_len; r++)
                    ChangeCoor(xyz_list1[r], RotMatix, TranVect, 
                        super_xyz_list1[r]);

                /* RMSD */
                rmsd=calRMSD(super_xyz_list1, xyz_list2);
                tmscore1=calTMscore(super_xyz_list1,xyz_list2,len1);
                tmscore2=calTMscore(super_xyz_list1,xyz_list2,len2);
                
                /* clean up */
                RotMatix.clear();
                TranVect.clear();
                super_xyz_list1.clear();
                xyz_list1.clear();
                xyz_list2.clear();
            }
            else if (aln_len!=0 && super_type>=2) // TM-score
            {
                if (input_mode>=6)
                    TMalign(aln1,aln2,pdb_entry1.chains[q],
                        pdb_entry1.chains[s],tmscore1,tmscore2,
                        (super_type==3),true);
                else
                    TMalign(aln1,aln2,pdb_entry1.chains[q],
                        pdb_entry2.chains[s],tmscore1,tmscore2,
                        (super_type==3),true);
            }

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
                if (super_type>=1)
                {
                    cout<<setiosflags(ios::fixed)<<setprecision(4)
                        <<"# IDali="<<1.*iden_len/aln_len
                        <<"\tidentity1="<<1.*iden_len/len1
                        <<"\tidentity2="<<1.*iden_len/len2
                        <<endl
                        <<"# Lali="<<aln_len<<"\tcoverage1="<<1.*aln_len/len1
                        <<"\tcoverage2="<<1.*aln_len/len2<<endl
                        <<"# RMSD="<<rmsd
                        <<"\tTM-score1="<<tmscore1
                        <<"\tTM-score2="<<tmscore2<<endl;
                }
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

                if (super_type>=1)
                {
                    cout<<"Aligned length= "<<aln_len
                        <<setiosflags(ios::fixed)<<setprecision(4)
                        <<", RMSD= "<<rmsd
                        <<", Seq_ID=n_identical/n_aligned= "
                        <<1.*iden_len/aln_len
                        <<"\nTM-score= "<<tmscore1
                        <<" (if normalized by length of Chain_1)"
                        <<"\nTM-score= "<<tmscore2
                        <<" (if normalized by length of Chain_2)\n"
                        <<endl;
                }

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
