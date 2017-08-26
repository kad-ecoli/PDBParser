/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef BackboneTorsion_HPP
#define BackboneTorsion_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

/* backbone torsion angles (Kappa, Alpha).
 * following DSSP's defination, kappa is supplementary angle
 * for the angle between i-2,i,i+2 CA atoms 
 *
 * if dist_instead == 1: calculate distance between
 * d(r-2,r), d(r-2,r+1), d(r-2,r+2), d(r-1,r+1), d(r-1,r+2), d(r,r+2) */
vector<vector<double> > KappaAlpha(ChainUnit& chain,int dist_instead=0)
{
    int L=chain.residues.size();
    vector<vector<double> > angle_mat;
    if (dist_instead==0)
    {
        vector<double> tmp_array(2,360.); // default torsion angles: 360
        angle_mat.assign(L,tmp_array);
    }
    else
    {
        vector<double> tmp_array(6,-1); // default distance: -1
        angle_mat.assign(L,tmp_array);
    }

    // coordinates of previous-previous residue
    vector<double> prev_prev_CA(3,0.);bool has_prev_prev_CA=false;
    // coordinates of previous-previous residue
    vector<double> prev_CA(3,0.);bool has_prev_CA=false;
    // coordinates of current residue
    vector<double> cur_CA(3,0.); bool has_cur_CA=false;
    // coordinates of next residue
    vector <double> next_CA(3,0); bool has_next_CA=false;
    // coordinates of next-next residue
    vector <double> next_next_CA(3,0); bool has_next_next_CA=false;

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

        if (dist_instead==0) // calculate kappa & alpha angles
        {
            if (has_prev_prev_CA && has_cur_CA && has_next_next_CA)
                angle_mat[r][0]=180.-rad2deg(Points2Angle(
                    prev_prev_CA,cur_CA,next_next_CA)); // kappa

            if (has_prev_CA && has_cur_CA && has_next_CA && has_next_next_CA)
                angle_mat[r][1]=rad2deg(Points2Dihedral(
                    prev_CA,cur_CA,next_CA,next_next_CA)); // alpha
        }
        else
        {
            if (has_prev_prev_CA && has_cur_CA)
                angle_mat[r][0]=Points2Distance(prev_prev_CA,cur_CA);
            if (has_prev_prev_CA && has_next_CA)
                angle_mat[r][1]=Points2Distance(prev_prev_CA,next_CA);
            if (has_prev_prev_CA && has_next_next_CA)
                angle_mat[r][2]=Points2Distance(prev_prev_CA,next_next_CA);
            if (has_prev_CA && has_next_CA)
                angle_mat[r][3]=Points2Distance(prev_CA,next_CA);
            if (has_prev_CA && has_next_next_CA)
                angle_mat[r][4]=Points2Distance(prev_CA,next_next_CA);
            if (has_cur_CA && has_next_next_CA)
                angle_mat[r][5]=Points2Distance(cur_CA,next_next_CA);
        }
    }
    
    prev_prev_CA.clear();
    prev_CA.clear();
    cur_CA.clear();
    next_CA.clear();
    next_next_CA.clear();
    return angle_mat;
}

/* backbone torsion angles (Omega, Phi, Psi) */
vector<vector<double> > OmegaPhiPsi(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: omega=phi=psi=360
    vector<double> tmp_array(3,360.);
    vector<vector<double> > angle_mat(L,tmp_array);

    // coordinates of previous residue
    vector<double> prev_CA(3,0.);bool has_prev_CA=false;
    vector<double> prev_C(3,0.); bool has_prev_C =false;
    // coordinates of current residue
    vector<double> cur_N(3,0.);  bool has_cur_N=false;
    vector<double> cur_CA(3,0.); bool has_cur_CA=false;
    vector<double> cur_C(3,0.);  bool has_cur_C=false;
    // coordinates of next residue
    vector <double> next_N(3,0); bool has_next_N=false;

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

/* secondary structure assignment by CA-CA distance 
 * 1 - coil, 2 - helix, 3 -turn, 4 - strand*/
int sec_str(const vector<double> dis_vec)
{
    for (int j=0;j<6;j++) if (dis_vec[j]<0) return 1; // coil
    double dis13=dis_vec[0];
    double dis14=dis_vec[1];
    double dis15=dis_vec[2];
    double dis24=dis_vec[3];
    double dis25=dis_vec[4];
    double dis35=dis_vec[5];

    double delta=2.1;
    if (fabs(dis15-6.37)<delta && fabs(dis14-5.18)<delta &&
        fabs(dis25-5.18)<delta && fabs(dis13-5.45)<delta &&
        fabs(dis24-5.45)<delta && fabs(dis35-5.45)<delta)
        return 2; //helix

    delta=1.42;
    if (fabs(dis15-13  )<delta && fabs(dis14-10.4)<delta &&
        fabs(dis25-10.4)<delta && fabs(dis13- 6.1)<delta &&
        fabs(dis24- 6.1)<delta && fabs(dis35- 6.1)<delta)
        return 4; //strand

    if(dis15 < 8)
        return 3; //turn

    return 1; // coil
}

#endif
