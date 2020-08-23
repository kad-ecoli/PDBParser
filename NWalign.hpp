/* header for Needleman-Wunsch global sequence alignment */
#ifndef NWalign_HPP
#define NWalign_HPP 1

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <ctype.h>

#include "PDBParser.hpp"
#include "FilePathParser.hpp"
#include "MathTools.hpp"
#include "ROTSUMalign.hpp"
#include "SarstAlign.hpp"
#include "SSalign.hpp"
#include "ThreeDblastAlign.hpp"

using namespace std;

const int gapopen_blosum62=-11;
const int gapext_blosum62=-1;

const int BLOSUM62[24][24]={
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},//A
{-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},//R
{-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},//N
{-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},//D
{ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},//C
{-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},//Q
{-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},//E
{ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},//G
{-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},//H
{-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},//I
{-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},//L
{-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},//K
{-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},//M
{-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},//F
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},//P
{ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},//S
{ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},//T
{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4},//W
{-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},//Y
{ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},//V
{-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},//B
{-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},//Z
{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},//X
{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1},//*
};


const int gapopen_blastn=-15;
const int gapext_blastn=-4;

const int blastn_matrix[24][24]={
//A  T  C  G  U
{ 2,-3,-3,-3,-3},//A
{-3, 2,-3,-3, 2},//T
{-3,-3, 2,-3,-3},//C
{-3,-3,-3, 2,-3},//G
{-3, 2,-3,-3, 2},//U
};

const string aa_list="ARNDCQEGHILKMFPSTWYVBZX*";
const string na_list="atcgu";

/* convert amino acid to int */
inline int aa2int(char aa)
{
    for (int i=0;i<aa_list.length();i++) if (aa_list[i]==aa) return i;
    if (aa!=toupper(aa)) return aa2int(toupper(aa));
    return aa_list.length();
}

inline int rna2int(char aa)
{
    for (int i=0;i<na_list.length();i++) if (na_list[i]==aa) return i;
    if (aa!=tolower(aa)) return aa2int(tolower(aa));
    return na_list.length();
}

vector<int> aa2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(aa2int(sequence[r]));
    return seq2int;
}

vector<int> rna2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(rna2int(sequence[r]));
    return seq2int;
}

/* read multiple-chain PDB and extract the sequence 
 * seq_type - 0: amino acid, 1: chi-1 rotamer sequence 
 *            2: SARST code, 3: 3d-blast sequence */
int read_pdb_as_fasta(const char *filename,vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list,
    ModelUnit &pdb_entry, const int seq_type=0)
{
    int atomic_detail=0; // only read CA
    if (seq_type==1 || seq_type==2) atomic_detail=2; // full atom structure
    int allowX=1;        // only allow ATOM and MSE

    string PDBid=basename_no_ext(filename);
    pdb_entry=read_pdb_structure(filename,atomic_detail,allowX);

    int seq_num=pdb_entry.chains.size();
    string sequence;
    vector<int> seq2int;
    int c,r;
    int empty_seq_num=0;
    for (c=0;c<seq_num;c++)
    {
        switch (seq_type)
        {
            case 1:sequence=getRotSeq(pdb_entry.chains[c]);break;
            case 2:sequence=pdb2sarst(pdb_entry.chains[c]);break;
            case 3:sequence=pdb2ThreeDblast(pdb_entry.chains[c]);break;
            case 4:sequence=pdb2ss(pdb_entry.chains[c]);break;
            default:sequence=pdb2fasta(pdb_entry.chains[c]);break;
        }

        if (sequence.length()==0)
        {
            empty_seq_num++;
            continue;
        }
        seq_list.push_back(sequence);

        switch (seq_type)
        {
            case 1:seq2int_list.push_back(RotSeq2int(sequence));break;
            case 2:seq2int_list.push_back(sarst2int(sequence));break;
            case 3:seq2int_list.push_back(ThreeDblast2int(sequence));break;
            case 4:seq2int_list.push_back(ss2int(sequence));break;
            case 5:seq2int_list.push_back(rna2int(sequence));break;
            default:seq2int_list.push_back(aa2int(sequence));break;
        }
        sequence.clear();
        seq2int.clear();
        name_list.push_back(PDBid+':'+pdb_entry.chains[c].chainID_full);
    }
    //pdb_entry.chains.clear();
    return (seq_num-empty_seq_num);
}

/* parse sequence typed by keyboard
 * seq_type - 0: amino acid, 1: chi-1 rotamer sequence
 *            2: SARST code, 3: 3d-blast sequence */
int get_stdin_seq(const char *stdin_seq, vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list,
    const int seq_type=0)
{
    string line=(string) stdin_seq;
    if (line.length()==0) return 0;
    string name=line.substr(0,10);

    vector<int> seq2int;
    string sequence;
    for (int i=0;i<line.length();i++)
    {
        sequence+=line[i];
        if (seq_type==1) // chi-1 rotamer sequence
            seq2int.push_back(RotSeq2int(line[i]));
        else if (seq_type==2) // SARST code
            seq2int.push_back(sarst2int(line[i]));
        else if (seq_type==3) // 3d-blast sequence
            seq2int.push_back(ThreeDblast2int(line[i]));
        else if (seq_type==4) // secondary structure
            seq2int.push_back(ss2int(line[i]));
        else // amino acid sequence
            seq2int.push_back(aa2int(line[i]));
    }

    name_list.push_back(name);
    seq_list.push_back(sequence);
    seq2int_list.push_back(seq2int);
    return 1;
}

/* read multiple-sequence fasta
 * seq_type - 0: amino acid, 1: chi-1 rotamer sequence
 *            2: SARST code, 3: 3d-blast sequence */
int read_fasta(const char *filename, vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list,
    const int seq_type=0)
{
    int seq_num=0;
    ifstream fp(filename, ios::in);
    int use_stdin=(strcmp(filename,"-")==0);
    if (!fp && !use_stdin)
    {
        cerr<<"ERROR! Cannot read file "<<filename<<endl;
        return 0;
    }

    string line,name,sequence;
    vector<int> seq2int; // aa2int
    int i;
    while (use_stdin?cin.good():fp.good())
    {
        use_stdin?getline(cin,line):getline(fp,line);
        if (line.empty()) continue;
       
        if (line[0]=='>')
        {
            seq_num++;
            if (seq_num>1)
            {
                seq_list.push_back(sequence);
                seq2int_list.push_back(seq2int);
                sequence.clear();
                seq2int.clear();
            }
            name.clear();
            for (i=1;i<line.length();i++)
            {
                if (isspace(line[i])) break;
                name+=line[i];
            }
            name_list.push_back(name);
        }
        else
        {
            for (i=0;i<line.length();i++)
            {
                sequence+=line[i];
                if (seq_type==1) // chi-1 rotamer sequence
                    seq2int.push_back(RotSeq2int(line[i]));
                else if (seq_type==2) // SARST code
                    seq2int.push_back(sarst2int(line[i]));
                else if (seq_type==3) // 3d-blast sequence
                    seq2int.push_back(ThreeDblast2int(line[i]));
                else if (seq_type==4) // secondary structure
                    seq2int.push_back(ss2int(line[i]));
                else if (seq_type==5) // RNA
                    seq2int.push_back(rna2int(line[i]));
                else // amino acid sequence
                    seq2int.push_back(aa2int(line[i]));
            }
        }
    }
    if (seq_num)
    {
        seq_list.push_back(sequence);
        seq2int_list.push_back(seq2int);
    }
    if (!use_stdin) fp.close();
    return seq_num;
}

/* initialize matrix in gotoh algorithm */
void init_gotoh_mat(vector<vector<int> >&JumpH, vector<vector<int> >&JumpV,
    vector<vector<int> >& P,vector<vector<int> >& S, 
    vector<vector<int> >& H, vector<vector<int> >& V,
    const int len1, const int len2, const int gapopen,const int gapext,
    const int glocal=0, const int alt_init=1)
{
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    for (i=0;i<len1+1;i++)
    {
        if (glocal<2) P[i][0]=4; // -
        JumpV[i][0]=i;
    }
    for (j=0;j<len2+1;j++)
    {
        if (glocal<1) P[0][j]=2; // |
        JumpH[0][j]=j;
    }
    if (glocal<2) for (i=1;i<len1+1;i++) S[i][0]=gapopen+gapext*(i-1);
    if (glocal<1) for (j=1;j<len2+1;j++) S[0][j]=gapopen+gapext*(j-1);
    if (alt_init==0)
    {
        for (i=1;i<len1+1;i++) H[i][0]=gapopen+gapext*(i-1);
        for (j=1;j<len2+1;j++) V[0][j]=gapopen+gapext*(j-1);
    }
    else
    {
        if (glocal<2) for (i=1;i<len1+1;i++) V[i][0]=gapopen+gapext*(i-1);
        if (glocal<1) for (j=1;j<len2+1;j++) H[0][j]=gapopen+gapext*(j-1);
        for (i=0;i<len1+1;i++) H[i][0]=-99999; // INT_MIN cause bug on ubuntu
        for (j=0;j<len2+1;j++) V[0][j]=-99999; // INT_MIN;
    }
}

/* locate the cell with highest alignment score. reset path after
 * the cell to zero */
void find_highest_align_score(
    const vector<vector<int> >& S, vector<vector<int> >& P,
    int &aln_score, const int len1,const int len2)
{
    // locate the cell with highest alignment score
    int max_aln_i=len1;
    int max_aln_j=len2;
    int i,j;
    for (i=0;i<len1+1;i++)
    {
        for (j=0;j<len2+1;j++)
        {
            if (S[i][j]>=aln_score)
            {
                max_aln_i=i;
                max_aln_j=j;
                aln_score=S[i][j];
            }
        }
    }

    // reset all path after [max_aln_i][max_aln_j]
    for (i=max_aln_i+1;i<len1+1;i++) for (j=0;j<len2+1;j++) P[i][j]=0;
    for (i=0;i<len1+1;i++) for (j=max_aln_j+1;j<len2+1;j++) P[i][j]=0;
}

/* calculate dynamic programming matrix using gotoh algorithm
 * S     - cumulative scorefor each cell
 * P     - string representation for path
 *         0 :   uninitialized, for gaps at N- & C- termini when glocal>0
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seq1)+1]*[len(seq2)+1]
 *
 * global - global or local alignment
 *         0 : global alignment (Needleman-Wunsch dynamic programming)
 *         1 : glocal-query alignment
 *         2 : glocal-both alignment
 *         3 : local alignment (Smith-Waterman dynamic programming)
 *
 * alt_init - whether to adopt alternative matrix initialization
 *         1 : use wei zheng's matrix initialization
 *         0 : use yang zhang's matrix initialization, does NOT work
 *             for glocal alignment
 */
int calculate_score_gotoh(
    const vector<int>& seq2int1, const vector<int>& seq2int2,
    vector<vector<int> >& JumpH, vector<vector<int> >& JumpV,
    vector<vector<int> >& P,const int ScoringMatrix[24][24],
    const int gapopen,const int gapext,const int glocal=0,
    const int alt_init=1)
{
    int len1=seq2int1.size();
    int len2=seq2int2.size();

    vector<int> temp_int(len2+1,0);
    vector<vector<int> > S(len1+1,temp_int);
    // penalty score for horizontal long gap
    vector<vector<int> > H(len1+1,temp_int);
    // penalty score for vertical long gap
    vector<vector<int> > V(len1+1,temp_int);
    
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    init_gotoh_mat(JumpH, JumpV, P, S, H, V, len1, len2,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score,left_score,up_score;
    for (i=1;i<len1+1;i++)
    {
        for (j=1;j<len2+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<len1 || glocal>=3)
            {
                H[i][j]=MAX(S[i][j-1]+gapopen,H[i][j-1]+gapext);
                JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            }
            else
            {
                H[i][j]=MAX(S[i][j-1],H[i][j-1]);
                JumpH[i][j]=(H[i][j]==H[i][j-1])?(JumpH[i][j-1]+1):1;
            }
            // penalty of consective insertion
            if (glocal<2 || j<len2 || glocal>=3)
            {
                V[i][j]=MAX(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=MAX(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]; // match-mismatch '\'
            if (seq2int1[i-1]<24 && seq2int2[j-1]<24)
                diag_score+=ScoringMatrix[seq2int1[i-1]][seq2int2[j-1]];
            left_score=H[i][j];     // deletion       '-'
            up_score  =V[i][j];     // insertion      '|'

            if (diag_score>=left_score && diag_score>=up_score)
            {
                S[i][j]=diag_score;
                P[i][j]+=1;
            }
            if (up_score>=diag_score && up_score>=left_score)
            {
                S[i][j]=up_score;
                P[i][j]+=2;
            }
            if (left_score>=diag_score && left_score>=up_score)
            {
                S[i][j]=left_score;
                P[i][j]+=4;
            }
            if (glocal>=3 && S[i][j]<0)
            {
                S[i][j]=0;
                P[i][j]=0;
                H[i][j]=0;
                V[i][j]=0;
                JumpH[i][j]=0;
                JumpV[i][j]=0;
            }
        }
    }
    int aln_score=S[len1][len2];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<len1+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<len2+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(S,P,aln_score,len1,len2);

    // release memory
    S.clear();
    H.clear();
    V.clear();
    return aln_score; // final alignment score
}

/* trace back dynamic programming path to diciper pairwise alignment */
void trace_back_gotoh(string seq1, string seq2,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& aln1, string& aln2)
{
    int len1=seq1.length();
    int len2=seq2.length();
    
    int i=len1;
    int j=len2;
    int gaplen,p;

    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            for (p=0;p<gaplen;p++) aln1='-'+aln1;
            aln2=seq2.substr(seq2.length()-gaplen,gaplen)+aln2;
            seq2=seq2.substr(0,seq2.length()-gaplen);
            j-=gaplen;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            aln1=seq1.substr(seq1.length()-gaplen,gaplen)+aln1;
            for (p=0;p<gaplen;p++) aln2='-'+aln2;
            seq1=seq1.substr(0,seq1.length()-gaplen);
            i-=gaplen;
        }
        else
        {
            if (i==0 && j!=0) // only in glocal alignment
            {
                aln2=seq2+aln2;
                for (p=0;p<seq2.length();p++) aln1='-'+aln1;
                break;
            }
            if (i!=0 && j==0) // only in glocal alignment
            {
                aln1=seq1+aln1;
                for (p=0;p<seq1.length();p++) aln2='-'+aln2;
                break;
            }
            aln1=seq1[seq1.length()-1]+aln1;
            aln2=seq2[seq2.length()-1]+aln2;
            seq1=seq1.substr(0,seq1.length()-1);
            seq2=seq2.substr(0,seq2.length()-1);
            i--;
            j--;
        }
    }   
}


/* trace back Smith-Waterman dynamic programming path to diciper 
 * pairwise local alignment */
void trace_back_sw(string seq1, string seq2,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& aln1, string& aln2)
{
    int len1=seq1.length();
    int len2=seq2.length();
    
    int i=len1;
    int j=len2;
    // find the first non-zero cell in P
    bool found_start_cell=false;
    for (i=len1;i>=0;i--)
    {
        for (j=len2;j>=0;j--)
        {
            if (P[i][j]!=0)
            {
                found_start_cell=true;
                break;
            }
        }
        if (found_start_cell) break;
    }
    if (i<0||j<0) return;
    seq1=seq1.substr(0,i);
    seq2=seq2.substr(0,j);

    int gaplen,p;
    while(P[i][j]!=0)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            for (p=0;p<gaplen;p++) aln1='-'+aln1;
            aln2=seq2.substr(seq2.length()-gaplen,gaplen)+aln2;
            seq2=seq2.substr(0,seq2.length()-gaplen);
            j-=gaplen;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            aln1=seq1.substr(seq1.length()-gaplen,gaplen)+aln1;
            for (p=0;p<gaplen;p++) aln2='-'+aln2;
            seq1=seq1.substr(0,seq1.length()-gaplen);
            i-=gaplen;
        }
        else
        {
            aln1=seq1[seq1.length()-1]+aln1;
            aln2=seq2[seq2.length()-1]+aln2;
            seq1=seq1.substr(0,seq1.length()-1);
            seq2=seq2.substr(0,seq2.length()-1);
            i--;
            j--;
        }
    }
}

/* entry function for NWalign 
 * seq_type - 0: amino acid, 1: chi1 rotamer, 2: sarst code, 
 *            3: 3d-blast sequence, 4 - ss */
int NWalign(const string& seq1, const string& seq2, 
    const vector<int>& seq2int1, const vector<int>& seq2int2, // aa2int
    string & aln1,string & aln2,const int seq_type=0, const int glocal=0)
{
    int len1=seq2int1.size();
    int len2=seq2int2.size();
    vector<int> temp_int(len2+1,0);
    vector<vector<int> > JumpH(len1+1,temp_int);
    vector<vector<int> > JumpV(len1+1,temp_int);
    vector<vector<int> > P(len1+1,temp_int);
    
    int aln_score;
    switch (seq_type)
    {
        case 1: // chi-1 rotamer
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
                ROTSUM8,gapopen_rotsum8,gapext_rotsum8,glocal);
            break;
        case 2: // sarst code
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
                BLOSUM62_sarst,gapopen_sarst,gapext_sarst,glocal);
            break;
        case 3: // 3d-blast sequence
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
                BLOSUM62_3dblast,gapopen_3dblast,gapext_3dblast,glocal);
            break;
        case 4: // ss, for some reason, using BLOSUM62 gives better result
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
                BLOSUM62_ss,gapopen_ss,gapext_ss,glocal);
            break;
        case 5: // rna
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,blastn_matrix,
                (glocal==3)?-5:gapopen_blastn, (glocal==3)?-2:gapext_blastn, glocal);
            break;
        default: // amino acid
            aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
                BLOSUM62,gapopen_blosum62,gapext_blosum62,glocal);
    }

    if (glocal<3)
        trace_back_gotoh(seq1,seq2,JumpH,JumpV,P,aln1,aln2);
    else
        trace_back_sw(seq1,seq2,JumpH,JumpV,P,aln1,aln2);

    JumpH.clear();
    JumpV.clear();
    P.clear();
    return aln_score; // aligment score
}

int get_seqID(const string& aln1, const string& aln2,
    string &aln_str,string &pos_str, int &iden_len,int &aln_len)
{
    iden_len=0;
    aln_len=0;
    for (int i=0;i<aln1.length();i++)
    {
        pos_str+=int('0')+(i+1)%10;
        if (aln1[i]==aln2[i] && aln1[i]!='-')
        {
            iden_len++;
            aln_str+=':';
        }
        else
        {
            aln_str+=' ';
        }
        aln_len+=(aln1[i]!='-' && aln2[i]!='-');
    }
    return iden_len;
}

#endif
