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

#define MAX(A,B) ((A)>(B)?(A):(B))

using namespace std;

const int gapopen_blosum62=-11;
const int gapext_blosum62=-1;

const int BLOSUM62[24][24]={
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
{ 4,-1,-2,-2 ,0,-1,-1 ,0,-2,-1,-1,-1,-1,-2,-1 ,1 ,0,-3,-2 ,0,-2,-1 ,0,-4},//A
{-1 ,5 ,0,-2,-3 ,1 ,0,-2 ,0,-3,-2 ,2,-1,-3,-2,-1,-1,-3,-2,-3,-1 ,0,-1,-4},//R
{-2 ,0 ,6 ,1,-3 ,0 ,0 ,0 ,1,-3,-3 ,0,-2,-3,-2 ,1 ,0,-4,-2,-3 ,3 ,0,-1,-4},//N
{-2,-2 ,1 ,6,-3 ,0 ,2,-1,-1,-3,-4,-1,-3,-3,-1 ,0,-1,-4,-3,-3 ,4 ,1,-1,-4},//D
{ 0,-3,-3,-3 ,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},//C
{-1 ,1 ,0 ,0,-3 ,5 ,2,-2 ,0,-3,-2 ,1 ,0,-3,-1 ,0,-1,-2,-1,-2 ,0 ,3,-1,-4},//Q
{-1 ,0 ,0 ,2,-4 ,2 ,5,-2 ,0,-3,-3 ,1,-2,-3,-1 ,0,-1,-3,-2,-2 ,1 ,4,-1,-4},//E
{ 0,-2 ,0,-1,-3,-2,-2 ,6,-2,-4,-4,-2,-3,-3,-2 ,0,-2,-2,-3,-3,-1,-2,-1,-4},//G
{-2 ,0 ,1,-1,-3 ,0 ,0,-2 ,8,-3,-3,-1,-2,-1,-2,-1,-2,-2 ,2,-3 ,0 ,0,-1,-4},//H
{-1,-3,-3,-3,-1,-3,-3,-4,-3 ,4 ,2,-3 ,1 ,0,-3,-2,-1,-3,-1 ,3,-3,-3,-1,-4},//I
{-1,-2,-3,-4,-1,-2,-3,-4,-3 ,2 ,4,-2 ,2 ,0,-3,-2,-1,-2,-1 ,1,-4,-3,-1,-4},//L
{-1 ,2 ,0,-1,-3 ,1 ,1,-2,-1,-3,-2 ,5,-1,-3,-1 ,0,-1,-3,-2,-2 ,0 ,1,-1,-4},//K
{-1,-1,-2,-3,-1 ,0,-2,-3,-2 ,1 ,2,-1 ,5 ,0,-2,-1,-1,-1,-1 ,1,-3,-1,-1,-4},//M
{-2,-3,-3,-3,-2,-3,-3,-3,-1 ,0 ,0,-3 ,0 ,6,-4,-2,-2 ,1 ,3,-1,-3,-3,-1,-4},//F
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4 ,7,-1,-1,-4,-3,-2,-2,-1,-2,-4},//P
{ 1,-1 ,1 ,0,-1 ,0 ,0 ,0,-1,-2,-2 ,0,-1,-2,-1 ,4 ,1,-3,-2,-2 ,0 ,0 ,0,-4},//S
{ 0,-1 ,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1 ,1 ,5,-2,-2 ,0,-1,-1 ,0,-4},//T
{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1 ,1,-4,-3,-2,11 ,2,-3,-4,-3,-2,-4},//W
{-2,-2,-2,-3,-2,-1,-2,-3 ,2,-1,-1,-2,-1 ,3,-3,-2,-2 ,2 ,7,-1,-3,-2,-1,-4},//Y
{ 0,-3,-3,-3,-1,-2,-2,-3,-3 ,3 ,1,-2 ,1,-1,-2,-2 ,0,-3,-1 ,4,-3,-2,-1,-4},//V
{-2,-1 ,3 ,4,-3 ,0 ,1,-1 ,0,-3,-4 ,0,-3,-3,-2 ,0,-1,-4,-3,-3 ,4 ,1,-1,-4},//B
{-1 ,0 ,0 ,1,-3 ,3 ,4,-2 ,0,-3,-3 ,1,-1,-3,-1 ,0,-1,-3,-2,-2 ,1 ,4,-1,-4},//Z
{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2 ,0 ,0,-2,-1,-1,-1,-1,-1,-4},//X
{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4 ,1},//*
};

const string aa_list="ARNDCQEGHILKMFPSTWYVBZX*";

/* convert amino acid to int */
inline int aa2int(char aa)
{
    for (int i=0;i<aa_list.length();i++) if (aa_list[i]==aa) return i;
    if (aa!=toupper(aa)) return aa2int(toupper(aa));
    return aa_list.length();
}

vector<int> aa2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(aa2int(sequence[r]));
    return seq2int;
}

/* read multiple-chain PDB and extract the sequence */
int read_pdb_as_fasta(const char *filename,vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list)
{
    int atomic_detail=0; // only read CA
    int allowX=1;        // only allow ATOM and MSE

    string PDBid=basename_no_ext(filename);
    ModelUnit pdb_entry=read_pdb_structure(filename,atomic_detail,allowX);

    int seq_num=pdb_entry.chains.size();
    string sequence;
    vector<int> seq2int;
    int c,r;
    for (c=0;c<seq_num;c++)
    {
        sequence=pdb2fasta(pdb_entry.chains[c]);
        seq_list.push_back(sequence);
        seq2int_list.push_back(aa2int(sequence));
        sequence.clear();
        seq2int.clear();
        name_list.push_back(PDBid+':'+pdb_entry.chains[c].chainID_full);
    }
    pdb_entry.chains.clear();
    return seq_num;
}

/* parse sequence typed by keyboard */
int get_stdin_seq(const char *stdin_seq, vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list)
{
    string line=(string) stdin_seq;
    if (line.length()==0) return 0;
    string name=line.substr(0,10);

    vector<int> seq2int;
    string sequence;
    for (int i=0;i<line.length();i++)
    {
         if (aa_list.find_first_of(toupper(line[i]))!=string::npos)
         {
             sequence+=line[i];
             seq2int.push_back(aa2int(line[i]));
         }
    }

    name_list.push_back(name);
    seq_list.push_back(sequence);
    seq2int_list.push_back(seq2int);
    return 1;
}

/* read multiple-sequence fasta */
int read_fasta(const char *filename, vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list)
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
                if (aa_list.find_first_of(toupper(line[i]))!=string::npos)
                {
                    sequence+=line[i];
                    seq2int.push_back(aa2int(line[i]));
                }
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

/* calculate dynamic programming matrix using gotoh algorithm
 * S     - cumulative scorefor each cell
 * P     - string representation for path 
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seq1)+1]*[len(seq2)+1]
 */
int calculate_score_gotoh(
    const vector<int>& seq2int1, const vector<int>& seq2int2,
    vector<vector<int> >& JumpH, vector<vector<int> >& JumpV,
    vector<vector<int> >& P,const int ScoringMatrix[24][24],
    const int gapopen,const int gapext)
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
    for (i=0;i<len1+1;i++)
    {
        P[i][0]=4;
        JumpV[i][0]=i;
    }
    for (j=0;j<len2+1;j++)
    {
        P[0][j]=2;
        JumpH[0][j]=j;
    }
    for (i=1;i<len1+1;i++)
    {
        S[i][0]=gapopen+gapext*(i-1);
        H[i][0]=gapopen+gapext*(i-1);
    }
    for (j=1;j<len2+1;j++)
    {
        S[0][j]=gapopen+gapext*(j-1);
        V[0][j]=gapopen+gapext*(j-1);
    }

    // fill S and P
    int diag_score,left_score,up_score;
    for (i=1;i<len1+1;i++)
    {
        for (j=1;j<len2+1;j++)
        {
            // penalty of consective deletion
            H[i][j]=MAX(S[i][j-1]+gapopen,H[i][j-1]+gapext);
            JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            // penalty of consective insertion
            V[i][j]=MAX(S[i-1][j]+gapopen,V[i-1][j]+gapext);
            JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;

            diag_score=S[i-1][j-1]+((seq2int1[i-1]<24&&seq2int2[j-1]<24)?
                ScoringMatrix[seq2int1[i-1]][seq2int2[j-1]]:0); // match '\'
            left_score=H[i][j];            // deletion       '-'
            up_score  =V[i][j];            // insertion      '|'

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
        }
    }
    int aln_score=S[len1][len2];
    S.clear();
    H.clear();
    V.clear();

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<len1;i++) P[i][0]=2; // |
    for (j=1;j<len2;j++) P[0][j]=4; // -
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
            aln1=seq1[seq1.length()-1]+aln1;
            aln2=seq2[seq2.length()-1]+aln2;
            seq1=seq1.substr(0,seq1.length()-1);
            seq2=seq2.substr(0,seq2.length()-1);
            i--;
            j--;
        }
    }   
}

/* entry function for NWalign */
int NWalign(const string& seq1, const string& seq2, 
    const vector<int>& seq2int1, const vector<int>& seq2int2, // aa2int
    string & aln1,string & aln2,const int ScoringMatrix[24][24],
    const int gapopen,const int gapext)
{
    int len1=seq2int1.size();
    int len2=seq2int2.size();
    vector<int> temp_int(len2+1,0);
    vector<vector<int> > JumpH(len1+1,temp_int);
    vector<vector<int> > JumpV(len1+1,temp_int);
    vector<vector<int> > P(len1+1,temp_int);

    int aln_score=calculate_score_gotoh(seq2int1,seq2int2,JumpH,JumpV,P,
        ScoringMatrix,gapopen,gapext);

    trace_back_gotoh(seq1,seq2,JumpH,JumpV,P,aln1,aln2);

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