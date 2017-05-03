/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

/* backbone torsion angles (Kappa, Alpha).
 * following DSSP's defination, kappa is supplementary angle
 * for the angle between i-2,i,i+2 CA atoms */
vector<vector<float> > KappaAlpha(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: kappa=alpha=360
    vector<float> tmp_array(2,360.);
    vector<vector<float> > angle_mat(L,tmp_array);

    // coordinates of previous-previous residue
    vector<float> prev_prev_CA(3,0.);bool has_prev_prev_CA=false;
    // coordinates of previous-previous residue
    vector<float> prev_CA(3,0.);bool has_prev_CA=false;
    // coordinates of current residue
    vector<float> cur_CA(3,0.); bool has_cur_CA=false;
    // coordinates of next residue
    vector <float> next_CA(3,0); bool has_next_CA=false;
    // coordinates of next-next residue
    vector <float> next_next_CA(3,0); bool has_next_next_CA=false;

    int r,a; // residue index, atom index
    
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_prev_prev_CA=false; has_prev_CA=false;
        has_cur_CA=false; has_next_CA=false; has_next_next_CA=false;

        if (r-2>=0) // find previous-previous residue atoms
        {
            for (a=0;a<chain.residues[r-2].atoms.size();a++)
            {
                if (chain.residues[r-2].atoms[a].name==" CA ")
                {
                    has_prev_prev_CA=true;
                    prev_prev_CA=chain.residues[r-2].atoms[a].xyz;
                }
            }
        }

        if (r-1>=0) // find previous residue atoms
        {
            for (a=0;a<chain.residues[r-1].atoms.size();a++)
            {
                if (chain.residues[r-1].atoms[a].name==" CA ")
                {
                    has_prev_CA=true;
                    prev_CA=chain.residues[r-1].atoms[a].xyz;
                }
            }
        }

        // find current residue atoms
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name==" CA ")
            {
                has_cur_CA=true;
                cur_CA=chain.residues[r].atoms[a].xyz;
            }
        }

        if (r+1<L) // find next residue atoms
        {
            for (a=0;a<chain.residues[r+1].atoms.size();a++)
            {
                if (chain.residues[r+1].atoms[a].name==" CA ")
                {
                    has_next_CA=true;
                    next_CA=chain.residues[r+1].atoms[a].xyz;
                }
            }
        }

        if (r+2<L) // find next-next residue atoms
        {
            for (a=0;a<chain.residues[r+2].atoms.size();a++)
            {
                if (chain.residues[r+2].atoms[a].name==" CA ")
                {
                    has_next_next_CA=true;
                    next_next_CA=chain.residues[r+2].atoms[a].xyz;
                }
            }
        }
 
        if (has_prev_prev_CA && has_cur_CA && has_next_next_CA)
            angle_mat[r][0]=180.-rad2deg(Points2Angle(
                prev_prev_CA,cur_CA,next_next_CA)); // kappa

        if (has_prev_CA && has_cur_CA && has_next_CA && has_next_next_CA)
            angle_mat[r][1]=rad2deg(Points2Dihedral(
                prev_CA,cur_CA,next_CA,next_next_CA)); // alpha
    }
    
    prev_prev_CA.clear();
    prev_CA.clear();
    cur_CA.clear();
    next_CA.clear();
    next_next_CA.clear();
    return angle_mat;
}

/* backbone torsion angles (Omega, Phi, Psi) */
vector<vector<float> > OmegaPhiPsi(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: omega=phi=psi=360
    vector<float> tmp_array(3,360.);
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

        if (r+1<L) // find next residue atoms
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
                prev_CA,prev_C,cur_N,cur_CA)); // omega

        if (has_prev_C && has_cur_N && has_cur_CA && has_cur_C)
            angle_mat[r][1]=rad2deg(Points2Dihedral(
                prev_C,cur_N,cur_CA,cur_C));   // phi

        if (has_cur_N && has_cur_CA && has_cur_C && has_next_N)
            angle_mat[r][2]=rad2deg(Points2Dihedral(
                cur_N,cur_CA,cur_C,next_N));   // psi
    }
    
    prev_CA.clear();
    prev_C.clear();
    cur_N.clear();
    cur_CA.clear();
    cur_C.clear();
    next_N.clear();
    return angle_mat;
}
