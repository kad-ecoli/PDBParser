/* header for sequence based structure alignment */
#include <string>
#include <vector>
#include <cmath>

#include "PDBParser.hpp"
#include "MathTools.hpp"
#include "GeometryTools.hpp"

using namespace std;

/* extract CA atom coordinates from chain1 and chain2 that are corresponding
 * to aligned regions shared by aln1 and aln2, return the aligned coordinates
 * xyz_list1 and xyz_list2
 * atomic_detail - 0 : assume first atom of a residue is CA
 *               >=1 : loop over all atoms of a residue to find CA
 */
void aln2coor(const string &aln1,const string &aln2,
    ChainUnit &chain1,ChainUnit &chain2,
    vector<vector<double> > &xyz_list1,vector<vector<double> > &xyz_list2,
    int atomic_detail=0)
{
    int L=MIN(aln1.length(),aln2.length());

    int pos; // pos is alignment position
    int r,r1=-1,r2=-1; // residue index of chain1 and chain2
    int a1=0,a2=0; // index for CA atom in current residue
    for (pos=0;pos<L;pos++)
    {
        if (aln1[pos]!='-')
        {
            for (r1+=1;r1<chain1.residues.size();r1++)
            {
                if (chain1.residues[r1].het==false)
                {
                    if (atomic_detail==0) break;
                    for (a1=0;a1<chain1.residues[r1].atoms.size();a1++)
                        if (chain1.residues[r1].atoms[a1].name==" CA ")
                            break;
                    if (a1<chain1.residues[r1].atoms.size())
                        break; // found 'ATOM' residue with CA atom
                }
            }
        }

        if (aln2[pos]!='-')
        {
            for (r2+=1;r2<chain2.residues.size();r2++)
            {
                if (chain2.residues[r2].het==false)
                {
                    if (atomic_detail==0) break;
                    for (a2=0;a2<chain2.residues[r2].atoms.size();a2++)
                        if (chain2.residues[r2].atoms[a2].name==" CA ")
                            break;
                    if (a2<chain2.residues[r2].atoms.size())
                        break; // found 'ATOM' residue with CA atom
                }
            }
        }

        if (aln1[pos]!='-' && aln2[pos]!='-')
        {
            xyz_list1.push_back(chain1.residues[r1].atoms[a1].xyz);
            xyz_list2.push_back(chain2.residues[r2].atoms[a2].xyz);
        }
    }
}

/* calculate RMSD of superposed structures */
double calRMSD(const vector<vector<double> >& xyz_list1, 
    const vector<vector<double> >& xyz_list2)
{
    double dist2=0; // sum of square difference
    int L=xyz_list1.size();
    for (int i=0;i<L;i++)
        dist2+=Points2Distance2(xyz_list1[i],xyz_list2[i]);

    double rms=0;
    if (L>0) rms=sqrt(dist2/L);
    return rms;
}


/* calculate TM-score of superposed structures, L is chain length */
double calTMscore(const vector<vector<double> >& xyz_list1, 
    const vector<vector<double> >& xyz_list2, int L)
{
    if (L<1) return 0.;

    double d0=(L>21)?(1.24*pow((L-15.),(1./3))-1.8):0.5;

    double d02=d0*d0;
    double TMscore=0; // sum of residue TMscore
    int Lali=xyz_list1.size();
    for (int i=0;i<Lali;i++)
        TMscore+=1./(1.+Points2Distance2(xyz_list1[i],xyz_list2[i])/d02);

    TMscore/=L;
    return TMscore;
}
