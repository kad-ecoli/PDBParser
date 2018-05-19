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
    /* one letter amino acid code */
    char  aa;

    /* residue non-specific angles and bonds */
    float CA_N_length ;
    float CA_C_length ;
    float C_O_length  ;
    float phi         ;
    float psi_im1     ;
    float omega       ;
    float peptide_bond;
    float CA_C_N_angle;
    float C_N_CA_angle;

    float CA_CB_length;
    float C_CA_CB_angle;

    /* residue specific angles */
    float N_CA_C_O_diangle;
    float CA_C_O_angle;
    float N_CA_C_angle;
    float N_C_CA_CB_diangle;

    /* initializer */
    Geo(const char aa='G')
    {
        CA_N_length      = 1.46;
        CA_C_length      = 1.52;
        C_O_length       = 1.23;
        phi              =-120.;
        psi_im1          = 140.;
        omega            =180.0;
        peptide_bond     = 1.33;
        CA_C_N_angle     =116.642992978143;
        C_N_CA_angle     =121.382215820277;
        CA_CB_length     = 1.52;
        C_CA_CB_angle    =109.5;

        N_CA_C_O_diangle =120.0; // LTRKDEQMHFYW
        CA_C_O_angle     =120.5;
        N_CA_C_angle     = 111.;
        N_C_CA_CB_diangle=122.6;
        assign(aa);
    }

    /* reassign to a new amino acid type */
    void assign(const char aa)
    {
        this->aa     =aa;
        if      (aa=='G')
        {
            N_CA_C_O_diangle=180.0;
            CA_C_O_angle=120.5117;
            N_CA_C_angle=110.8914;
        }
        else if (aa=='A')
        {
            N_CA_C_O_diangle=-60.5;
            N_CA_C_angle=111.068;
            N_C_CA_CB_diangle=122.6860;
        }
        else if (aa=='P')
        {
            N_CA_C_O_diangle=-45.0;
            CA_C_O_angle=120.2945;
            N_CA_C_angle=112.7499;
            N_C_CA_CB_diangle=115.2975;
        }
        else if (string("SCUVIN").find(aa)!=string::npos)
        {
            N_CA_C_O_diangle= -60.0;
            if (aa=='S')
            {
                N_CA_C_angle=111.2812;
                N_C_CA_CB_diangle=122.6618;
            }
            else if (aa=='C' || aa=='U')
            {
                N_CA_C_angle= 110.8856;
                N_C_CA_CB_diangle=122.5037;
            }
            else if (aa=='V')
            {
                CA_C_O_angle=120.5686;
                N_CA_C_angle=109.7698;
                N_C_CA_CB_diangle=123.2347;
            }
            else if (aa=='I')
            {
                CA_C_O_angle=120.5403;
                N_CA_C_angle=109.7202;
                N_C_CA_CB_diangle=123.2347;
            }
            else if (aa=='N') 
            {
                CA_C_O_angle=120.4826;
                N_CA_C_angle=111.5;
                N_C_CA_CB_diangle=123.2254;
            }
        }
        else if (aa=='L')
        {
            CA_C_O_angle=120.4647;
            N_CA_C_angle=110.8652;
            N_C_CA_CB_diangle=122.4948;
        }
        else if (aa=='T')
        {
            CA_C_O_angle=120.5359;
            N_CA_C_angle=110.7014;
            N_C_CA_CB_diangle=123.0953;
        }
        else if (aa=='R'||aa=='K'||aa=='O')
        {
            CA_C_O_angle=120.54;
            N_CA_C_angle=111.08;
            N_C_CA_CB_diangle=122.76;
            if (aa=='R') N_CA_C_angle=110.98;
        }
        else if (aa=='D')
        {
            CA_C_O_angle=120.51;
            N_CA_C_angle=111.03;
            N_C_CA_CB_diangle=122.82;
        }
        else if (aa=='Q')
        {
            CA_C_O_angle=120.5029;
            N_CA_C_angle=111.0849;
            N_C_CA_CB_diangle=122.8134;
        }
        else if (aa=='E')
        {
            CA_C_O_angle=120.511;
            N_CA_C_angle=111.1703;
            N_C_CA_CB_diangle=122.8702;
        }
        else if (aa=='M')
        {
            CA_C_O_angle=120.4816;
            N_CA_C_angle=110.9416;
            N_C_CA_CB_diangle=122.6733;
        }
        else if (aa=='H')
        {
            CA_C_O_angle=120.4732;
            N_CA_C_angle=111.0859;
            N_C_CA_CB_diangle=122.6711;
        }
        else if (aa=='F')
        {
            CA_C_O_angle=120.5316;
            N_CA_C_angle=110.7528;
            N_C_CA_CB_diangle=122.6054;
        }
        else if (aa=='Y')
        {
            CA_C_O_angle=120.5434;
            N_CA_C_angle=110.9288;
            N_C_CA_CB_diangle=122.6023;
        }
        else if (aa=='W')
        {
            CA_C_O_angle=120.5117;
            N_CA_C_angle=110.8914;
            N_C_CA_CB_diangle=122.6112;
        }
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

/* Adds a residue to given chain. "residue" is the last added
 * residue. "geo" must be pre-assigned according to the residue
 * to be added. */
void add_residue_from_geo(ChainUnit &chain, ResidueUnit &residue,
    const Geo &geo)
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
        else if (residue.atoms[a].name==" C  ")
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

    /* correct old Carbonyl oxygen atom */
    vector <float> prevO_coord(3,0);
    calculateCoordinates(prevO_coord, N_coord, prevCA_coord, 
        prevC_coord, geo.C_O_length, geo.CA_C_O_angle, 180.0);
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" O  ")
        {
            for (int i=0;i<3;i++)
                chain.residues.back().atoms[a].xyz[i]=prevO_coord[i];
        }
    }
    prevO_coord.clear();

    /* add new residue to chain */
    residue.resi++;
    makeRes(residue, N_coord, CA_coord, C_coord, O_coord, geo);
    chain.residues.push_back(residue);

    /* clean up */
    prevN_coord.clear();
    prevCA_coord.clear();
    prevC_coord.clear();
}

/* add OXT to last residue. "residue" must be the last residue.
 * "geo" is reassigned according to "residue" */
void add_terminal_oxygen(ChainUnit &chain, ResidueUnit &residue, 
    Geo &geo, const float psi)
{
    geo.assign(aa3to1(residue.resn));

    /* last residue */
    vector <float> N_coord;
    vector <float> CA_coord;
    vector <float> C_coord;
    int a;
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" N  ") 
            N_coord=residue.atoms[a].xyz;
        else if (residue.atoms[a].name==" CA ")
            CA_coord=residue.atoms[a].xyz;
        else if (residue.atoms[a].name==" C  ")
            C_coord=residue.atoms[a].xyz;
    }

    /* OXT coordinates */
    AtomUnit atom;
    atom.name=" OXT";
    atom.xyz.assign(3,0);
    calculateCoordinates(atom.xyz, N_coord, CA_coord, C_coord,
        geo.C_O_length, geo.CA_C_N_angle, psi);

    chain.residues.back().atoms.push_back(atom);

    /* clean up */
    N_coord.clear();
    CA_coord.clear();
    C_coord.clear();
    atom.name.clear();
    atom.xyz.clear();
}

/* add C beta given backbone N C CA */
void addCbeta(ResidueUnit &residue, Geo &geo)
{
    geo.assign(aa3to1(residue.resn));

    /* check if C beta already exist*/
    int a;
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" CB ")
            break; // do not add CB if there is already one
    }
    if (a<residue.atoms.size()) return; // there is already a C beta

    /* get backbone coordinates */
    vector<float> N_coord;
    vector<float> C_coord;
    vector<float> CA_coord;
    for (a=0;a<residue.atoms.size();a++)
    {
        if (residue.atoms[a].name==" N  ") 
            N_coord=residue.atoms[a].xyz;
        else if (residue.atoms[a].name==" CA ")
            CA_coord=residue.atoms[a].xyz;
        else if (residue.atoms[a].name==" C  ")
            C_coord=residue.atoms[a].xyz;
    }

    /* get CB coordinates */
    AtomUnit atom;
    atom.name=" CB ";
    atom.xyz.assign(3,0);
    calculateCoordinates(atom.xyz, N_coord, C_coord, CA_coord, 
        geo.CA_CB_length, geo.C_CA_CB_angle, geo.N_C_CA_CB_diangle);
    residue.atoms.push_back(atom);

    /* clean up */
    atom.xyz.clear();
    atom.name.clear();
    N_coord.clear();
    C_coord.clear();
    CA_coord.clear();
}

void addCbeta(ChainUnit &chain, Geo &geo, const bool gly_beta=false)
{
    for (int r=0;r<chain.residues.size();r++)
    {
        if (gly_beta || chain.residues[r].resn!="GLY")
            addCbeta(chain.residues[r], geo);
    }
}

void addCbeta(ModelUnit &model, Geo &geo, const bool gly_beta=false)
{
    for (int c=0;c<model.chains.size();c++)
        addCbeta(model.chains[c], geo, gly_beta);
}

/* add_cb - whether to add CB atoms. 
 *          0: backbone only
 *          1: backbone + CB for non gly residue
 *          2: backbone + CB for all residue, including gly
 */
void make_chain(ChainUnit &chain, const string sequence, 
    const vector<vector<float> >&rama_table, const bool add_oxt=true,
    const int add_cb=1)
{
    if (sequence.size()==0) return;

    chain.chainID='A';
    chain.chainID_full="A";

    Geo geo(sequence[0]);
    ResidueUnit residue;
    initialize_res(residue, geo);
    chain.residues.push_back(residue);

    int r;
    for (r=1;r<sequence.size();r++)
    {
        geo.assign(sequence[r]);

        geo.omega=rama_table[r][0];
        geo.phi=rama_table[r][1];
        geo.psi_im1=rama_table[r-1][2];

        add_residue_from_geo(chain, residue, geo);
    }

    if (add_cb) addCbeta(chain, geo, (add_cb==2));

    if (add_oxt)
        add_terminal_oxygen(chain, residue, geo, rama_table[r-1][2]);

    /* clean up */
    residue.resn.clear();
    residue.atoms.clear();
}
