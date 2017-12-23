/* ce2tm: convert ccealign superposition to TM-align alignment.
 * The original ccealign module prefers low RMSD over long coverage,
 * causing alignments to have short alignments. This module uses
 * dynamic programming to obtain TM-align style alignment so that
 * alignment length is longer.
 */
#ifndef ce2tm_HPP
#define ce2tm_HPP 1

#include "MathTools.hpp"
#include "GeometryTools.hpp"
#include "PDBParser.hpp"
#include "NWalign.hpp"

using namespace std;

const float gapopen_tm=-1;
const float gapext_tm=-1;

/* initialize matrix in gotoh algorithm */
void init_gotoh_mat(vector<vector<int> >&JumpH, vector<vector<int> >&JumpV,
    vector<vector<int> >& P,vector<vector<float> >& S, 
    vector<vector<float> >& H, vector<vector<float> >& V,
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
    const vector<vector<float> >& S, vector<vector<int> >& P,
    float &aln_score, const int len1,const int len2)
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
float calculate_score_gotoh(
    vector<vector<int> >& JumpH, vector<vector<int> >& JumpV,
    vector<vector<int> >& P,const vector<vector<float> > TM_matrix,
    const float gapopen,const float gapext,const int glocal=0,
    const int alt_init=1)
{
    int len1=TM_matrix.size();
    int len2=TM_matrix[0].size();

    vector<float> temp_float(len2+1,0);
    vector<vector<float> > S(len1+1,temp_float);
    // penalty score for horizontal long gap
    vector<vector<float> > H(len1+1,temp_float);
    // penalty score for vertical long gap
    vector<vector<float> > V(len1+1,temp_float);
    
    // fill first row/colum of JumpH,jumpV and path matrix P
    init_gotoh_mat(JumpH, JumpV, P, S, H, V, len1, len2,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int i,j;
    float diag_score,left_score,up_score;
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

            diag_score=S[i-1][j-1]+TM_matrix[i-1][j-1]; // aligned '\'
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
    float aln_score=S[len1][len2];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<len1+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<len2+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3) find_highest_align_score(S,P,aln_score,len1,len2);

    // release memory
    S.clear();
    H.clear();
    V.clear();
    return aln_score; // final alignment score
}

/* entry function for TM-score based dynamic programming */
float NWalign(const string& seq1, const string& seq2, 
    const vector<vector<float> > TM_matrix,
    string & aln1,string & aln2, const int glocal=0)
{
    int len1=TM_matrix.size();
    int len2=TM_matrix[0].size();
    vector<int> temp_int(len2+1,0);
    vector<vector<int> > JumpH(len1+1,temp_int);
    vector<vector<int> > JumpV(len1+1,temp_int);
    vector<vector<int> > P(len1+1,temp_int);
    
    float aln_score=calculate_score_gotoh(JumpH,JumpV,P,TM_matrix,
        gapopen_tm,gapext_tm,glocal);

    if (glocal<3) trace_back_gotoh(seq1,seq2,JumpH,JumpV,P,aln1,aln2);
    else trace_back_sw(seq1,seq2,JumpH,JumpV,P,aln1,aln2);

    JumpH.clear();
    JumpV.clear();
    P.clear();
    return aln_score; // aligment score
}

/* Transform xyz into super_xyz using transform matrix rVal.
 * The PyMOL-specific transform matrix rVal consists of
 * [1] a 3x3 matrix containing the rotation in the upper-left quadrant;
 * [2] a 1x3 translation to be applied before rotation in the bottom row 
 *     (matrix[12],matrix[13],matrix[14]).
 * [3] a 3x1 translation to be applied after rotation in the right-hand
 *     column (matrix[3],matrix[7],matrix[11])
 *
 * So, you can translate+rotate+translate with this one matrix.
 * In other words, if the matrix is:
 * [ m0  m1  m2  m3
 *   m4  m5  m6  m7
 *   m8  m9 m10 m11
 *   m12 m13 m14 m15 ]
 *
 * xyz=[x0,x1,x2] will be transformed into super_xyz=[y0,y1,y2] as follows:
 * y0 = m0*(x0+m12) + m1*(x1+m13) +  m2*(x2+m14) + m3
 * y1 = m4*(x0+m12) + m5*(x1+m13) +  m6*(x2+m14) + m7
 * y2 = m8*(x0+m12) + m9*(x1+m13) + m10*(x2+m14) + m11 
 */
void transform_selection(vector<float> xyz, 
    const vector<vector<float> > &rVal, vector<float> &super_xyz)
{
    xyz[0]+=rVal[3][0];
    xyz[1]+=rVal[3][1];
    xyz[2]+=rVal[3][2];
    super_xyz[0]=rVal[0][0]*xyz[0]+rVal[0][1]*xyz[1]+rVal[0][2]*xyz[2]+rVal[0][3];
    super_xyz[1]=rVal[1][0]*xyz[0]+rVal[1][1]*xyz[1]+rVal[1][2]*xyz[2]+rVal[1][3];
    super_xyz[2]=rVal[2][0]*xyz[0]+rVal[2][1]*xyz[1]+rVal[2][2]*xyz[2]+rVal[2][3];
}

/* transform the whole L*3 coordinate matrix */
void transform_selection(vector<vector<float> >xyz,
    const vector<vector<float> > &rVal, vector<vector<float> >&super_xyz)
{
    for (int r=0;r<xyz.size();r++)
        transform_selection(xyz[r], rVal, super_xyz[r]);
}

/* calculate algnment from superposed coordinates 
 * aln_moving & aln_fixed should be unaligned sequences
 */
void super2aln(const vector<vector<float> >&coord_moving,
    const vector<vector<float> >&coord_fixed,
    string &aln_moving, string &aln_fixed, string &aln_str, float &rmsd, 
    float &tmscore_moving, float &tmscore_fixed, int &L_ali, float &seqID)
{
    int L_moving=coord_moving.size();
    int L_fixed=coord_fixed.size();
    int L=MIN(L_moving,L_fixed);

    float d0_moving=(L_moving>21)?(1.24*pow((L_moving-15.),(1./3))-1.8):0.5;
    float d0_fixed=(L_fixed>21)?(1.24*pow((L_fixed-15.),(1./3))-1.8):0.5;
    float d02_moving=d0_moving*d0_moving;
    float d02_fixed=d0_fixed*d0_fixed;
    float d02=MIN(d02_moving,d02_fixed);

    /* calculate TM-score matrix for each pair of residues */
    vector<float> tmp_array(L_fixed,0);
    vector<vector<float> > TM_matrix(L_moving,tmp_array);
    int i,j;
    for (i=0;i<L_moving;i++)
    {
        for (j=0;j<L_fixed;j++)
        {
            TM_matrix[i][j]=1./(1.+Points2Distance2(
                coord_moving[i],coord_fixed[j])/d02);
        }
    }

    /* unaligned sequence */
    string seq_moving=aln_moving;
    string seq_fixed=aln_fixed;
    aln_moving.clear();
    aln_fixed.clear();

    /* make alignment */
    NWalign(seq_moving, seq_fixed, TM_matrix, aln_moving, aln_fixed);

    /* calculate TM-score */
    rmsd=0;
    L_ali=0;
    tmscore_moving=0;
    tmscore_fixed=0;
    seqID=0;
    int moving_idx=-1;
    int fixed_idx=-1;
    float di2;
    aln_str.clear();
    for (i=0;i<aln_moving.size();i++)
    {
        moving_idx+=(aln_moving[i]!='-');
        fixed_idx+=(aln_fixed[i]!='-');
        if (aln_moving[i]=='-' || aln_fixed[i]=='-') 
        {
            aln_str+=' ';
            continue;
        }

        L_ali+=1;
        di2=Points2Distance2(coord_moving[moving_idx],coord_fixed[fixed_idx]);
        rmsd+=di2;
        tmscore_moving+=1./(1.+di2/d02_moving);
        tmscore_fixed+=1./(1.+di2/d02_fixed);
        seqID+=1.*(aln_moving[i]==aln_fixed[i]);
        aln_str+=(di2<5*5)?':':'.';
    }
    rmsd=sqrt(rmsd/L_ali);
    seqID/=L_ali;
    tmscore_moving/=L_moving;
    tmscore_fixed/=L_fixed;
}

/* calculate aignment from superposition matrix and unsuperposed chains */
void super2aln(ChainUnit &chain_moving, ChainUnit &chain_fixed,
    const vector<vector<float> >&rVal, string &aln_moving, string &aln_fixed,
    string &aln_str, float &rmsd, float &tmscore_moving, float &tmscore_fixed, 
    int &L_ali, float &seqID)
{
    vector<vector<float> > coord_moving=getCoords(chain_moving);
    vector<vector<float> > coord_fixed=getCoords(chain_fixed);
    vector<vector<float> > coord_super=coord_moving;

    transform_selection(coord_moving, rVal, coord_super);
    if (chain_moving.sequence.length()==0) pdb2fasta(chain_moving);
    if (chain_fixed.sequence.length()==0)  pdb2fasta(chain_fixed);
    aln_moving=chain_moving.sequence;
    aln_fixed=chain_fixed.sequence;

    super2aln(coord_super, coord_fixed, aln_moving, aln_fixed, aln_str,
        rmsd, tmscore_moving, tmscore_fixed, L_ali, seqID);
}

#endif
