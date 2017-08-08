/* convert backbone torsion angles to structural alphabet */
#ifndef StructuralAphabet_HPP
#define StructuralAphabet_HPP 1
#include <cstring>
#include "BackboneTorsion.hpp"

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

const char* ThreeDblast_matrix=
    "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
    "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
    "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
    "ZZZZZZZZZZZZZZZZZZLZZZZZZZZZZZZZZZZZ"
    "ZZZZZZZZSSSSSWWLWLLLLLIILZZZZZZZZZZZ"
    "ZZZZZZSSSSSSWWWWWWLLLLLIILLLRZZZZZZZ"
    "ZZZZZSSSSSSSSWWWWWWWLDACIILLQRRZZZZZ"
    "ZZZSSSSSSSSSSWWWVVVVMDDBDLQQQQQRZZZZ"
    "QSSSSSSSSSSVVVVVVVVVVMGGGGQQQQQQQQZZ"
    "PPSSSSTTTTTTTVVVVVVVVMGMGGQQQQQQQQRQ"
    "PPPPPPTTTTTTVVVVVVVMMMMMMQQQQQQRRQPP"
    "PPPPPTTTTTTTTTVVVVVXMXMMMMXXXRRRRRPP"
    "PPPPTTTTNNNTXXXXXXXXXXXXXXXXXXXXXRRP"
    "NPTTTNNKKKKKKKXKXXXXZXXXXXXXXXXXXNNN"
    "HNNNKKKKKKKKKZZZZZZZZZZZZXZXXXXXHHHH"
    "HHHKFFKKKKKZZZZZZZZZZZZZZZZZXXXXHHHH"
    "EEEFFFKKZZZZZZZZZZZZZZZZZZZZZXXXXNHH"
    "EEFFFFZZZZZZZZZZZZZZZZZZZZZZZZZZNNHH";

const char* ss_vec="XCHTE"; // undefined, coil, helix, strand, turn

/* convert pdb chain to secondary structure */
string pdb2ss(ChainUnit& chain)
{
    vector<vector<float> > dist_mat=KappaAlpha(chain,1);
    string ss_seq="";
    int L=chain.residues.size();
    for (int r=0;r<L;r++)
        ss_seq+=ss_vec[sec_str(dist_mat[r])];
    return ss_seq;
}

/* convert pdb entry to secondary structure
 * ShowSeqLen - whether to show residue number for each chain */
string pdb2ss(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=1)
{
    stringstream buf;
    string ss_seq="";
    for (int c=0;c<pep.chains.size();c++)
    {
        ss_seq=pdb2ss(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID_full;
        if (ShowSeqLen) buf<<'\t'<<ss_seq.length();
        buf<<'\n'<<ss_seq<<'\n';
    }
    return buf.str();
}

/* convert pdb chain to SARST (Structural similarity search Aided 
 * by Ramachandran Sequential Transformation) code */
string pdb2sarst(ChainUnit& chain)
{
    vector<vector<float> > angle_mat=OmegaPhiPsi(chain);
    chain.sarst="";
    int L=chain.residues.size();
    int dim=36; //int(sqrt(strlen(sarst_matrix)));

    int i,j; // coordinate in sarst_matrix
    for (int r=0;r<L;r++)
    {
        if (has_atom_name(chain.residues[r])==0) continue;
        if (angle_mat[r][1]==360 || angle_mat[r][2]==360)
        {
            chain.sarst+='X';
        }
        else
        {
            i=int((180-angle_mat[r][2])/(360/dim)); // psi
            j=int((180+angle_mat[r][1])/(360/dim)); // phi
            chain.sarst+=sarst_matrix[i*dim+j];
        }
    }
    angle_mat.clear();
    return chain.sarst;
}

/* convert pdb entry to SARST (Structural similarity search Aided 
 * by Ramachandran Sequential Transformation) code
 * ShowSeqLen - whether to show residue number for each chain */
string pdb2sarst(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=1)
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
    return buf.str();
}

/* convert pdb chain to structural alphabet in 3D-blast */
string pdb2ThreeDblast(ChainUnit& chain)
{
    vector<vector<float> > angle_mat=KappaAlpha(chain);
    string ThreeDblast_seq="";
    int L=chain.residues.size();
    int dim=36; //int(sqrt(strlen(sarst_matrix)));

    int i,j; // coordinate in sarst_matrix
    char SA_code;
    float kappa,alpha;
    for (int r=0;r<L;r++)
    {
        if (has_atom_name(chain.residues[r])==0) continue;
        kappa=angle_mat[r][0];
        alpha=angle_mat[r][1];
        if (kappa==360 || alpha==360)
        {
            ThreeDblast_seq+='U'; // chain termini
        }
        else
        {
            i=int((180-kappa)/(360/dim)); // kappa
            j=int((180+alpha)/(360/dim)); // alpha
            SA_code=ThreeDblast_matrix[i*dim+j];

            if (SA_code=='A' && kappa<=114 && alpha>=46 )
                ThreeDblast_seq+='Y'; // unpublished detail
            else
                ThreeDblast_seq+=SA_code;
        }
    }
    angle_mat.clear();
    return ThreeDblast_seq;
}

/* convert pdb entry to structure alphabet in 3D-blast 
 * ShowSeqLen - whether to show residue number for each chain */
string pdb2ThreeDblast(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=1)
{
    stringstream buf;
    string ThreeDblast_seq="";
    for (int c=0;c<pep.chains.size();c++)
    {
        ThreeDblast_seq=pdb2ThreeDblast(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID_full;
        if (ShowSeqLen) buf<<'\t'<<ThreeDblast_seq.length();
        buf<<'\n'<<ThreeDblast_seq<<'\n';
    }
    return buf.str();
}

#endif
