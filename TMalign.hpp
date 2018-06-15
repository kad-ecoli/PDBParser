#include "PDBParser.hpp"
#include "TMalign.h"

using namespace std;

/* adaptor function to convert ChainUnit into TMalign style data structure */
int read_PDB(const ChainUnit& chain, double **xa, char *seq, int *resno)
{
    int r,a; // residue, atom
    for (r=0;r<chain.residues.size();r++)
    {
        seq[r]=aa3to1(chain.residues[r].resn);
        //resno[r]=chain.residues[r].resi;
        resno[r]=r;

        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name==" CA ")
            {
                xa[r][0]=chain.residues[r].atoms[a].xyz[0];
                xa[r][1]=chain.residues[r].atoms[a].xyz[1];
                xa[r][2]=chain.residues[r].atoms[a].xyz[2];
                break;
            }
        }
    }
    seq[r]='\0';
    return r;
}

/* free alignment */
void TMalign(string & aln1, string & aln2,
    ChainUnit & chain1, const ChainUnit & chain2,
    float & tmscore1, float & tmscore2, bool fast_opt, bool I_opt)
{
    tmscore1=tmscore2=0;

    /* extract alignment if necessary */
    vector<string> sequence; // get value from alignment file
    if (I_opt)
    {
        sequence.push_back(aln1);
        sequence.push_back(aln2);
    }
    aln1.clear();
    aln2.clear();

    /* dummy variables */
    bool i_opt = false; // flag for -i, with user given initial alignment
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    double Lnorm_ass=0;
    double d0_scale =0;

    /* variables for storing chain data */
    int xlen=chain1.residues.size();
    int ylen=chain2.residues.size();
    char   *seqx, *seqy;     // for the protein sequence 
    int    *secx, *secy;     // for the secondary structure 
    int    *xresno, *yresno; // residue number for fragment gapless threading
    double **xa, **ya;       // for input vectors xa[0...xlen-1][0..2] and
                             // ya[0...ylen-1][0..2], in general,
                             // ya is regarded as native structure 
                             // --> superpose xa onto ya
    
    /* variables for chain 1 */
    NewArray(&xa, xlen, 3);
    seqx = new char[xlen + 1];
    secx = new int[xlen];
    xresno = new int[xlen];
    read_PDB(chain1, xa, seqx, xresno);
    make_sec(xa, xlen, secx); // secondary structure assignment

    /* variables for chain 2 */
    NewArray(&ya, ylen, 3);
    seqy = new char[ylen + 1];
    yresno = new int[ylen];
    secy = new int[ylen];
    read_PDB(chain2, ya, seqy, yresno);
    make_sec(ya, ylen, secy);

    /* variables for a pair of alignment */
    double t0[3], u0[3][3];
    double TM1, TM2;
    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
    double d0_0, TM_0;
    double d0A, d0B, d0u, d0a;
    double d0_out=5.0;
    string seqM;              // for output alignment
    double rmsd0 = 0.0;
    int L_ali;                // Aligned length in standard_TMscore
    double Liden=0;
    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
    int n_ali=0;
    int n_ali8=0;

    /* entry function for structure alignment */
    TMalign_main(
        xa, ya, xresno, yresno, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
        seqM, aln1, aln2,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        i_opt, I_opt, a_opt, u_opt, d_opt, fast_opt);

    /* return value */
    tmscore2=TM1;
    tmscore1=TM2;

    /* clear memory */
    seqM.clear();
    DeleteArray(&ya, ylen);
    delete [] seqy;
    delete [] secy;
    delete [] yresno;
    DeleteArray(&xa, xlen);
    delete [] seqx;
    delete [] secx;
    delete [] xresno;
    return;
}
