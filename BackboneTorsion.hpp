/* calculate backbone torsion angles (omega, psi, phi) */
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

const char* sarst_matrix=
    "SSSIIIIQQQQQQQQZZZZZZZRRYYYYYYYYSSSS"
    "SSIIIIIQQQQFFFFZZZZZZZZYYYYYYYYYYSSS"
    "SIIIIIIILFFFFFFFZZZZZZZYYYYYYYYYYSSS"
    "SSIIIGGHHFFFFFFFZZZZZZZZYYYYYYYYYSSS"
    "SMIGGGGHHLLFFLFLZZZZZZZYYYYYYYYYZZZS"
    "MMGMGGHHHLLLLLLLLZZZZZZZYZYYYYYZZYSM"
    "MMMMMMMHHLLLLLLLZZZZZZZYYYZYYYZZYZZM"
    "MMMMMMMHHLLLLLLZZZZZZZZPYYYZZZZZZZMM"
    "MMMMMMMMLLLLLLLZZZZZZPPPPPYYYZYZZZZM"
    "MMMMMMMMMLLLLLZZZZZZZPPPPPZZZZZZZZZM"
    "MMMMMMMMMMLLLZZZZPZZPPPPPPZPZZZZZZZZ"
    "MMMMMMMMMNNNZZZZZZZZPPPPPPPPZZZZZZZM"
    "MMMMNNNNNNNNNZZZZZZZPPPPPPPZZPZWWWZM"
    "NNNNNNNNNNNKZZZZZZZPPPPPPPPPPPZZZZZW"
    "NNNNNNNNNNKNZZZZZZZPPPPPPPPPPWWZZZZZ"
    "WNNNNNNNNNKKZZZZZZZPPPPPPPPPWWWWZZZZ"
    "ZNNNNNNNNKKKKZZZZZZZPPPPPPPWWWWWWZZZ"
    "ZNNNNNNNNKKKKKZZZZZZZPPPPPPWWWWWWWZZ"
    "ZNNNNNNNKKKKDDDZZZZZZZPPPPWWWWWWWWWZ"
    "ZNNNNNNNKKKKDDDDZZZZZZZZPWWWWWWWWWWZ"
    "WVVNNNNNKKEDDDDDDZZZZZZZZWWWWWWWWWWZ"
    "VVVVVVNEEEBBDDDDDZZZZZZZWWWWWWWWWWWW"
    "VVVVVVVEEEEACCCCDZZZZZZZWWWWWWWWWZWV"
    "VVVVVVVEEEECCCTTTZZZZZZZWWWZWWWWZZWV"
    "VVVVVVVVEEEETTTTTTZZZZZRRWWZZZZZZZZV"
    "VVVVVVVVEEETTTTTTTZZZZZZRRZZZZZZZZZV"
    "VVVVVVVVVTTTTTTTTTZZZZZRRRRZZZZZZZZZ"
    "ZVVVVVVVVTTTTTTTTZZZZZRRRZRZZZZZZZZV"
    "VZVVVVVVVTTTTTTTZZZZZRRRRRZZZRZZZZZZ"
    "ZZZVVVVVVTTTTTTZZZZZRRRRRRRZZZZZZZZS"
    "ZZSVVVVQQQZZZZZZZZZRRRRRRRRRZZZZZZZZ"
    "SSZSSVQQQQQQZZZZZZZRRRRRRRRRZZZZZSSS"
    "SSSSSQQQQQQQQZZZZZZRRRRRRRRRRZYSSSSS"
    "SSSSSQQQQQQQQQZZZZZZRRRRRRRRYYYSSSSS"
    "SSSSSQQQQQQQQQQZZZZZZRRRRRRYYYYYSSSS"
    "SSSSIIQQQQQQQQQZZZZZZZRRRYYYYYYYSSSS";

/* BackboneTorsion - backbone torsion angles (Omega, Phi, Psi) */
vector<vector<float> > BackboneTorsion(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: omega=180, phi=psi=360
    vector<float> tmp_array(3,360.); tmp_array[0]=180.;
    vector<vector<float> > angle_mat(L,tmp_array);

    // coordinates of previous residue
    vector<float> prev_CA(3,0.);bool has_prev_CA=false;
    vector<float> prev_C(3,0.); bool has_prev_C =false;
    // coordinates of current residue
    vector<float> cur_N(3,0.);  bool has_cur_N=false;
    vector<float> cur_CA(3,0.); bool has_cur_CA=false;
    vector<float> cur_C(3,0.);  bool has_cur_C=false;
    // coordinates of next residue
    vector <float> next_N(3,0); bool has_next_N=false;

    int r,a; // residue index, atom index
    
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_prev_CA=false; has_prev_C=false;
        has_cur_N=false;   has_cur_CA=false; has_cur_C=false;
        has_next_N=false;
        
        if (r>0) // find previous residue atoms
        {
            for (a=0;a<chain.residues[r-1].atoms.size();a++)
            {
                if (chain.residues[r-1].atoms[a].name==" CA ")
                {
                    has_prev_CA=true;
                    prev_CA=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C  ")
                {
                    has_prev_C=true;
                    prev_C=chain.residues[r-1].atoms[a].xyz;
                }
            }
        }

        // find current residue atoms
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name==" N  ")
            {
                has_cur_N=true;
                cur_N=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CA ")
            {
                has_cur_CA=true;
                cur_CA=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C  ")
            {
                has_cur_C=true;
                cur_C=chain.residues[r].atoms[a].xyz;
            }
        }

        // find next residue atoms
        if (r+1<L)
        {
            for (a=0;a<chain.residues[r+1].atoms.size();a++)
            {
                if (chain.residues[r+1].atoms[a].name==" N  ")
                {
                    has_next_N=true;
                    next_N=chain.residues[r+1].atoms[a].xyz;
                }
            }
        }
 
        if (has_prev_CA && has_prev_C && has_cur_N && has_cur_CA)
            angle_mat[r][0]=rad2deg(Points2Dihedral(
                prev_CA,prev_C,cur_N,cur_CA));

        if (has_prev_C && has_cur_N && has_cur_CA && has_cur_C)
            angle_mat[r][1]=rad2deg(Points2Dihedral(
                prev_C,cur_N,cur_CA,cur_C));

        if (has_cur_N && has_cur_CA && has_cur_C && has_next_N)
            angle_mat[r][2]=rad2deg(Points2Dihedral(
                cur_N,cur_CA,cur_C,next_N));
    }
    
    prev_CA.clear();
    prev_C.clear();
    cur_N.clear();
    cur_CA.clear();
    cur_C.clear();
    next_N.clear();
    return angle_mat;
}

/* convert pdb chain to SARST (Structural similarity search Aided 
 * by Ramachandran Sequential Transformation) code */
string pdb2sarst(ChainUnit& chain)
{
    vector<vector<float> > angle_mat=BackboneTorsion(chain);
    string sarst_seq;
    int L=chain.residues.size();
    int dim=int(sqrt(strlen(sarst_matrix)));

    int i,j; // coordinate in sarst_matrix
    for (int r=0;r<L;r++)
    {
        if (angle_mat[r][1]==360 || angle_mat[r][2]==360)
        {
            sarst_seq+='X';
        }
        else
        {
            i=int((180-angle_mat[r][2])/(360/dim));
            j=int((180+angle_mat[r][1])/(360/dim));
            sarst_seq+=sarst_matrix[i*dim+j];
        }
    }
    return sarst_seq;
}

/* convert pdb entry to SARST (Structural similarity search Aided 
 * by Ramachandran Sequential Transformation) code
 * ShowSeqLen - whether to show residue number for each chain */
string pdb2sarst(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=0)
{
    stringstream buf;
    string sarst_seq="";
    for (int c=0;c<pep.chains.size();c++)
    {
        sarst_seq=pdb2sarst(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID_full;
        if (ShowSeqLen) buf<<'\t'<<sarst_seq.length();
        buf<<'\n'<<sarst_seq<<'\n';
    }
    sarst_seq.clear();
    return buf.str();
}
