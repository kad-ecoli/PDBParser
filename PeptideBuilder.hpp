/* Re-implementation of the "PeptideBuilder" algorithm described by:
 * Tien, Matthew Z., et al. "PeptideBuilder: A simple Python library to
 * generate model peptides." PeerJ 1 (2013): e80. */
#include <string>
#include <vector>
#include <cmath>

#include "PDBParser.hpp"
#include "MathTools.hpp"
#include "GeometryTools.hpp"

using namespace std;

/* define backbone geometry of generic amino acid */
class Geo
{
    public:
    float CA_N_length ;
    float CA_C_length ;
    float N_CA_C_angle;
    float C_O_length  ;
    float CA_C_O_angle;
    float phi         ;
    float psi_im1     ;
    float omega       ;
    float peptide_bond;
    float CA_C_N_angle;
    float C_N_CA_angle;
    char  aa; // one letter amino acid code
    float N_CA_C_O_diangle;

    /* initializer */
    Geo(const char aa)
    {
        CA_N_length  = 1.46;
        CA_C_length  = 1.52;
        N_CA_C_angle = 111.;
        C_O_length   = 1.23;
        CA_C_O_angle =120.5;
        phi          =-120.;
        psi_im1      = 140.;
        omega        =180.0;
        peptide_bond = 1.33;
        CA_C_N_angle =116.642992978143;
        C_N_CA_angle =121.382215820277;
        assign(aa);
    }

    /* reassign to a new amino acid type */
    void assign(const char aa)
    {
        this->aa     =aa;
        N_CA_C_O_diangle =120.0; // LTRKDEQMHFYW
        if      (aa=='G') N_CA_C_O_diangle=180.0;
        else if (aa=='A') N_CA_C_O_diangle=-60.5;
        else if (aa=='P') N_CA_C_O_diangle=-45.0;
        else if (string("SCUVIN").find(aa)!=string::npos)
            N_CA_C_O_diangle= -60.0;
    }
};

void makeRes(ResidueUnit &residue, const vector<float> &N_coord, 
    const vector<float> &CA_coord, const vector<float> &C_coord,
    const vector<float> &O_coord, const Geo &geo)
{
    residue.resn=aa1to3(geo.aa);
    residue.atoms.clear();

    AtomUnit atom;
    atom.xyz.assign(3,0);
    int i=0;

    atom.name=" N  ";
    for (i=0;i<3;i++) atom.xyz[i]=N_coord[i];
    residue.atoms.push_back(atom);

    atom.name=" CA ";
    for (i=0;i<3;i++) atom.xyz[i]=CA_coord[i];
    residue.atoms.push_back(atom);

    atom.name=" C  ";
    for (i=0;i<3;i++) atom.xyz[i]=C_coord[i];
    residue.atoms.push_back(atom);

    atom.name=" O  ";
    for (i=0;i<3;i++) atom.xyz[i]=O_coord[i];
    residue.atoms.push_back(atom);

    /* clean up */
    atom.xyz.clear();
    atom.name.clear();
}

void calculateCoordinates(vector<float>&D, 
    const vector<float>&refA, const vector<float>&refB, const vector<float>&refC, 
    const float l, const float ang, const float di)
{
    vector<float>CA(3,0);
    vector<float>CB(3,0);
    subtract(refA,refB,CA);
    subtract(refB,refC,CB);

    /* Plane Parameters */
    float A=(CA[1]*CB[2]) - (CA[2]*CB[1]);
    float B=(CA[2]*CB[0]) - (CA[0]*CB[2]);
    float G=(CA[0]*CB[1]) - (CA[1]*CB[0]);

    /* Dot Product Constant */
    float F=sqrt(CB[0]*CB[0] + CB[1]*CB[1] + CB[2]*CB[2]
        ) * l * cos(deg2rad(ang));

    float cons=B*CB[2]-CB[1]*G;
    cons=sqrt(cons*cons*(-(F*F)*(A*A+B*B+G*G)+(B*B*(CB[0]*CB[0]+CB[2]*CB[2]
        ) + A*A*(CB[1]*CB[1]+CB[2]*CB[2])- (2*A*CB[0]*CB[2]*G) + (
        CB[0]*CB[0]+ CB[1]*CB[1])*G*G- (2*B*CB[1])*(A*CB[0]+CB[2]*G))*l*l));
    float denom=B*B*(CB[0]*CB[0]+CB[2]*CB[2]) + A*A*(CB[1]*CB[1]+CB[2]*CB[2]
        ) - (2*A*CB[0]*CB[2]*G) + G*G*(CB[0]*CB[0]+CB[1]*CB[1]
        ) - (2*B*CB[1])*(A*CB[0]+CB[2]*G);

    float X=((B*B*CB[0]*F)-(A*B*CB[1]*F)+(F*G)*(-A*CB[2]+CB[0]*G)+cons
        )/denom;

    float Y,Z;
    if ((B==0 or CB[2]==0) && (CB[1]==0 or G==0))
    {
        float const1=sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(l-X)*(l+X)));
        Y= ((-A*B*X)+const1)/(B*B+G*G);
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G));
    }
    else
    {
        Y= ((A*A*CB[1]*F)*(B*CB[2]-CB[1]*G)+ G*(-F*pow(B*CB[2]-CB[1]*G,2
            ) + CB[0]*cons) - A*( B*B*CB[0]*CB[2]*F- B*CB[0]*CB[1]*F*G + \
            CB[2]*cons)) / ((B*CB[2]-CB[1]*G)*denom);
        Z= ((A*A*CB[2]*F)*(B*CB[2]-CB[1]*G) + (B*F)*pow(B*CB[2]-CB[1]*G,2
            ) + (A*CB[0]*F*G)*(-B*CB[2]+CB[1]*G) - B*CB[0]*cons + \
            A*CB[1]*cons) / ((B*CB[2]-CB[1]*G)*denom);
    }

    
    /* GET THE NEW VECTOR from the orgin */
    vector<float>tmp_array(3,0);
    tmp_array[0]=X;
    tmp_array[1]=Y;
    tmp_array[2]=Z;
    vectorsum(tmp_array, refC, D);

    float angle=di-rad2deg(Points2Dihedral(refA, refB, refC, D));
    CoordinateRotation(D,refC,refB, angle, tmp_array);
    for (int i=0;i<3;i++) D[i]=tmp_array[i];

    /* clean up */
    tmp_array.clear();
    CA.clear();
    CB.clear();
}

void initialize_res(ResidueUnit &residue, const Geo &geo)
{
    vector<float> CA_coord(3,0);

    vector<float> C_coord(3,0);
    C_coord[0]=geo.CA_C_length;

    vector<float> N_coord(3,0);
    N_coord[0]=geo.CA_N_length * cos(deg2rad(geo.N_CA_C_angle));
    N_coord[1]=geo.CA_N_length * sin(deg2rad(geo.N_CA_C_angle));

    /* Create Carbonyl oxygen atom (to be moved later) */
    vector<float> O_coord(3,0);
    calculateCoordinates(O_coord, N_coord, CA_coord, C_coord, 
        geo.C_O_length, geo.CA_C_O_angle, geo.N_CA_C_O_diangle);

    residue.resi=1;
    residue.icode=' ';
    residue.het=false;
    makeRes(residue, N_coord, CA_coord, C_coord, O_coord, geo);

    /* clean up */
    CA_coord.clear();
    C_coord.clear();
    N_coord.clear();
    O_coord.clear();
}

/* Adds a residue to chain A model 0 of the given structure, and
 * returns the new structure. The residue to be added is determined by
 * the geometry object given as second argument.
 */
void add_residue_from_geo(ChainUnit &chain, ResidueUnit &residue, const Geo &geo)
{
    /* getReferenceResidue */
    vector <float> prevN_coord;
    vector <float> prevCA_coord;
    vector <float> prevC_coord;
    int a;
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" N  ") 
            prevN_coord=residue.atoms[a].xyz;
        else if (residue.atoms[a].name==" CA ")
            prevCA_coord=residue.atoms[a].xyz;
        if (residue.atoms[a].name==" C  ")
            prevC_coord=residue.atoms[a].xyz;
    }

    /* new coordinates */
    vector<float> CA_coord(3,0);
    vector<float> C_coord(3,0);
    vector<float> N_coord(3,0);

    calculateCoordinates(N_coord, prevN_coord, prevCA_coord, prevC_coord,
        geo.peptide_bond, geo.CA_C_N_angle, geo.psi_im1);
    calculateCoordinates(CA_coord, prevCA_coord, prevC_coord, N_coord,
        geo.CA_N_length, geo.C_N_CA_angle, geo.omega);
    calculateCoordinates(C_coord, prevC_coord, N_coord, CA_coord, 
        geo.CA_C_length, geo.N_CA_C_angle, geo.phi);

    /* Create Carbonyl oxygen atom (to be moved later) */
    vector<float> O_coord(3,0);
    calculateCoordinates(O_coord, N_coord, CA_coord, C_coord, 
        geo.C_O_length, geo.CA_C_O_angle, geo.N_CA_C_O_diangle);
    //calculateCoordinates(O_coord, N_coord, CA_coord, C_coord,
        //geo.C_O_length, geo.CA_C_O_angle, 180.0);

    /* correct old Carbonyl oxygen atom */
    vector <float> prevO_coord(3,0);
    calculateCoordinates(prevO_coord, N_coord, prevCA_coord, prevC_coord,
        geo.C_O_length, geo.CA_C_O_angle, 180.0);
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" O  ")
        {
            for (int i=0;i<3;i++)
            {
                chain.residues.back().atoms[a].xyz[i]=prevO_coord[i];
            }
        }
    }

    residue.resi++;
    makeRes(residue, N_coord, CA_coord, C_coord, O_coord, geo);
    chain.residues.push_back(residue);

    /* clean up */
    prevN_coord.clear();
    prevCA_coord.clear();
    prevC_coord.clear();
}

void make_chain(ChainUnit &chain, const string sequence, 
    const vector<vector<float> >&rama_table)
{
    if (sequence.size()==0) return;

    chain.chainID='A';
    chain.chainID_full="A";

    Geo geo(sequence[0]);
    ResidueUnit residue;
    initialize_res(residue, geo);
    chain.residues.push_back(residue);

    cout<<geo.aa<<endl;
    for (int r=1;r<sequence.size();r++)
    {
        geo.assign(sequence[r]);

        geo.omega=rama_table[r][0];
        geo.phi=rama_table[r][1];
        geo.psi_im1=rama_table[r-1][2];

        add_residue_from_geo(chain, residue, geo);
    }

    /* clean up */
    residue.resn.clear();
    residue.atoms.clear();
}
