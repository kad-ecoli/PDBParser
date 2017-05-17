/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

/* sidechain torsion angles (chi1, chi2, chi3, chi4) */
vector<vector<float> >SidechainTorsion(ChainUnit& chain)
{
    int L=chain.residues.size();
    // default torsion angles: 360
    vector<float> tmp_array(4,360.);
    vector<vector<float> >angle_mat(L,tmp_array);

    // for chi1
    vector<float> N(3,0.);  bool has_N=false;
    vector<float> CA(3,0.); bool has_CA=false;
    vector<float> CB(3,0.); bool has_CB=false;
    vector<float> SG(3,0.); bool has_SG=false; // CYS
    vector<float> CG1(3,0.);bool has_CG1=false;// ILE VAL
    //vector<float> CG2(3,0.);bool has_CG2=false;// VAL, fake chi2
    vector<float> OG(3,0.); bool has_OG=false; // SER
    vector<float> OG1(3,0.);bool has_OG1=false;// THR
    vector<float> CG(3,0.); bool has_CG=false; // ARG ASN ASP GLN GLU HIS LEU
                                               // LYS MET PHE PRO TRP TYR
    // for chi2
    vector<float> CD(3,0.);bool has_CD=false;  // ILE ARG LYS GLU GLN
    vector<float> CD1(3,0.);bool has_CD1=false;// LEU
    vector<float> SD(3,0.);bool has_SD=false;  // MET

    // for chi3
    vector<float> CE(3,0.);bool has_CE=false;  // MET LYS
    vector<float> NE(3,0.);bool has_NE=false;  // ARG

    // for chi4
    vector<float> CZ(3,0.);bool has_CZ=false;  // ARG
    vector<float> NZ(3,0.);bool has_NZ=false;  // LYS

    int r,a; // residue index, atom index
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_N=false;  has_CA=false; has_CB=false; has_CG=false; // chi1
        has_CG1=false;has_OG=false; has_OG1=false;// chi1

        has_CD=false; has_CD1=false;has_SD=false; // chi2
        //has_CG2=false; // fake val chi2

        has_CE=false; has_NE=false; // chi3

        has_CZ=false; has_NZ=false; // chi4

        // find current residue atoms
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            // chi1
            if (chain.residues[r].atoms[a].name==" N  ")
            {
                has_N=true;
                N=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CA ")
            {
                has_CA=true;
                CA=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CB ")
            {
                has_CB=true;
                CB=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CG ")
            {
                has_CG=true;
                CG=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" SG ")
            {
                has_SG=true;
                SG=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CG1")
            {
                has_CG1=true;
                CG1=chain.residues[r].atoms[a].xyz;
            }
            //else if (chain.residues[r].atoms[a].name==" CG2") // fake CG2
            //{
                //has_CG2=true;
                //CG2=chain.residues[r].atoms[a].xyz;
            //}
            else if (chain.residues[r].atoms[a].name==" OG ")
            {
                has_OG=true;
                OG=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" OG1")
            {
                has_OG1=true;
                OG1=chain.residues[r].atoms[a].xyz;
            }
            // chi2
            else if (chain.residues[r].atoms[a].name==" CD ")
            {
                has_CD=true;
                CD=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" CD1")
            {
                has_CD1=true;
                CD1=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" SD ")
            {
                has_SD=true;
                SD=chain.residues[r].atoms[a].xyz;
            }
            // chi3
            else if (chain.residues[r].atoms[a].name==" CE ")
            {
                has_CE=true;
                CE=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" NE ")
            {
                has_NE=true;
                NE=chain.residues[r].atoms[a].xyz;
            }
            // chi4
            else if (chain.residues[r].atoms[a].name==" CZ ")
            {
                has_CZ=true;
                CZ=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" NZ ")
            {
                has_NZ=true;
                NZ=chain.residues[r].atoms[a].xyz;
            }
        }

        // calculate chi1
        if (!(has_N && has_CA && has_CB)) continue;
        if (chain.residues[r].resn=="CYS")
        {
            if (has_SG) angle_mat[r][0]=rad2deg(Points2Dihedral(N,CA,CB,SG));
            continue;
        }
        else if ((chain.residues[r].resn=="ILE"||chain.residues[r].resn=="VAL")
                && has_CG1)
        {
            angle_mat[r][0]=rad2deg(Points2Dihedral(N,CA,CB,CG1));
            if (chain.residues[r].resn=="VAL")
            {
                //if (has_CG2) 
                    //angle_mat[r][1]=rad2deg(Points2Dihedral(N,CA,CB,CG2));
                continue;
            }
        }
        else if (chain.residues[r].resn=="SER")
        {
            if (has_OG) angle_mat[r][0]=rad2deg(Points2Dihedral(N,CA,CB,OG));
            continue;
        }
        else if (chain.residues[r].resn=="THR")
        {
            if (has_OG1) angle_mat[r][0]=rad2deg(Points2Dihedral(N,CA,CB,OG1));
            continue;
        }
        else
        {
            if (!has_CG) continue;
            angle_mat[r][0]=rad2deg(Points2Dihedral(N,CA,CB,CG));
            if (chain.residues[r].resn=="PRO") continue; // skip pro chi2
        }

        // calculate chi2
        if (chain.residues[r].resn=="LEU")
        {
            if (has_CD1) angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG,CD1));
            continue;
        }
        else if (chain.residues[r].resn=="MET")
        {
            if (has_SD) angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG,SD));
            else continue;
        }
        else if (chain.residues[r].resn=="ILE")
        {
            if (has_CG1 && has_CD) 
                angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG1,CD));
            else if (has_CG1 && has_CD1) 
                angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG1,CD1));
            else if (has_CG && has_CD1)
                angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG1,CD));
            continue;
        }
        else
        {
            if (!has_CD) continue;
            angle_mat[r][1]=rad2deg(Points2Dihedral(CA,CB,CG,CD));
        }

        // calculate chi3
        if (chain.residues[r].resn=="MET")
        {
            if (has_CE) angle_mat[r][2]=rad2deg(Points2Dihedral(CB,CG,SD,CE));
            continue;
        }
        if (chain.residues[r].resn=="ARG")
        {
            if (has_NE) angle_mat[r][2]=rad2deg(Points2Dihedral(CB,CG,CD,NE));
            else continue;
        }
        else 
        {
            if (!has_CE) continue;
            angle_mat[r][2]=rad2deg(Points2Dihedral(CB,CG,CD,CE));
        }

        if (chain.residues[r].resn=="ARG" && has_CZ) 
            angle_mat[r][3]=rad2deg(Points2Dihedral(CG,CD,NE,CZ));
        else if (chain.residues[r].resn=="LYS" && has_CE)
            angle_mat[r][3]=rad2deg(Points2Dihedral(CG,CD,CE,NZ));
    }
    
    // clean chi1
    N.clear(); CA.clear(); CB.clear();
    SG.clear();CG1.clear();OG.clear();OG1.clear();CG.clear();
    // clean chi2
    CD.clear();CD1.clear();SD.clear();
    // clean chi3
    CE.clear();NE.clear();
    // clean chi4
    CZ.clear();NZ.clear();
    return angle_mat;
}
