/* parse PDB file into data structure similar to Bio.PDB in biopython
 * (model - chain - residue - atom). */
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "pstream.h"
#include "GeometryTools.hpp"

using namespace std;

struct AtomUnit    // struct for each atom entry
{
    string name;       // atom name
    vector<float> xyz; // coordinate
};

struct ResidueUnit // struct for each residue
{
    bool het;               // true - HETATM, false - ATOM
    int resi;               // residue sequence number
    char icode;             // insertion code
    string resn;            // residue name
    vector<AtomUnit> atoms; // list of atoms
};

struct ChainUnit  // struct for each chain
{
    string chainID_full;          // chain ID, might be more than 1 char
    char chainID;                 // short chain ID, must be 1 char
    string sequence;              // sequence converted from CA coordinate
    vector<ResidueUnit> residues; // list of residues
};

struct ModelUnit  // struct for each model in mult-model PDB
{
    vector<ChainUnit> chains; // list of chains
};

/* parse one line in PDB file, append the data to pep. 
 * used by read_pdb_structure
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET
 *         2 - all, converting MSE to MET, 3 - all, no conversion
 *
 * return 0 if line not parsed, 1 if line is parsed
 */
int parse_pdb_line(const string line,ModelUnit &pep, ChainUnit &chain,
    ResidueUnit &residue, AtomUnit &atom, map<char,string> &chainIDmap,
    const int atomic_detail=2,const int allowX=1)
{
    string record_name=line.substr(0,6);
    char altLoc=line[16];

    atom.name=line.substr(12,4);
    residue.resn=line.substr(17,3);

    if ((allowX==0 && record_name!="ATOM  ")||
        (allowX==1 && record_name!="ATOM  " &&  
            !(record_name=="HETATM" && residue.resn=="MSE"))||
        (allowX>=2 && record_name!="ATOM  " && record_name!="HETATM"))
        return 0;
    
    // ignore alternatively locating residues
    if (altLoc!=' ' && altLoc!='A') return 0;

    if ((atomic_detail==0 && atom.name!=" CA ")||
        (atomic_detail==1 && atom.name!=" CA " && atom.name!=" N  " && 
         atom.name!=" C  " && atom.name!=" O  ")) return 0;
    
    if (residue.resn=="MSE" && allowX<3)
    {
        record_name="ATOM  ";
        residue.resn="MET";
    }
    residue.het=(record_name=="HETATM");
    chain.chainID=line[21];
    if (chainIDmap.find(chain.chainID)==chainIDmap.end())
        chain.chainID_full=chain.chainID;
    else
        chain.chainID_full=chainIDmap[chain.chainID];
    residue.resi=atoi(line.substr(22,4).c_str());
    residue.icode=line[26];
    atom.xyz[0]=atof(line.substr(30,8).c_str());
    atom.xyz[1]=atof(line.substr(38,8).c_str());
    atom.xyz[2]=atof(line.substr(46,8).c_str());

    int chain_index=-1;
    for (int c=0;c<pep.chains.size();c++)
        if (pep.chains[c].chainID_full==chain.chainID_full) chain_index=c;
    if (chain_index==-1)
    {
        pep.chains.push_back(chain);
        chain_index=pep.chains.size()-1;
    }

    if (pep.chains[chain_index].residues.size()==0||
        pep.chains[chain_index].residues.back().resi !=residue.resi||
        pep.chains[chain_index].residues.back().icode!=residue.icode)
        pep.chains[chain_index].residues.push_back(residue);
    
    pep.chains[chain_index].residues.back().atoms.push_back(atom);
    return 1;
}

/* atomic_detail: 0 - CA only, 1 - backbone heavy atoms (CA C N O), 2 - all atom
 * allowX: 0 - ATOM, 1 - ATOM and MSE, converting MSE to MET, 
 *         2 - all, converting MSE to MET, 3 - all, no conversion
 * filename: full filename path, stdin if filename=="-"
 */
ModelUnit read_pdb_structure(const char *filename,
    const int atomic_detail=2,const int allowX=1)
{
    ModelUnit pep;

    string line="";
    string record_name="ATOM  ";
    char altLoc=' ';

    AtomUnit atom;
    atom.xyz.assign(3,0);

    ResidueUnit residue;

    ChainUnit chain;

    map<char,string> chainIDmap;
    string filename_str=(string) filename;
    
    if (filename_str.length()>=18 &&
        filename_str.substr(filename_str.length()-18,18)=="-pdb-bundle.tar.gz")
    {
        // best effort/minimal tarball
        redi::ipstream fp_gz("tar -xOzf "+filename_str+
            " --wildcards *-chain-id-mapping.txt");
        map <string,map<char,string> > PDBmap;
        vector<string> PDBvec;

        string PDBfile=""; // PDB format file in tarball
        char chainID;
        string chainID_full;

        while(fp_gz.good())
        {
            getline(fp_gz,line);
            if (line.length()==0 || (PDBfile.length()==0 && line[0]==' '))
                continue;
            if (line[0]!=' ') // new PDBfile
            {
                PDBfile=line.substr(0,20);
                PDBmap[PDBfile]=chainIDmap;
                PDBvec.push_back(PDBfile);
            }
            else // new chain
            {
                istringstream ss(line);
                ss>>chainID>>chainID_full;
                PDBmap[PDBfile][chainID]=chainID_full;
            }
        }
        fp_gz.close();

        int i,c;
        for (i=0;i<PDBvec.size();i++)
        {
            PDBfile=PDBvec[i];
            redi::ipstream fp_gz2("tar -xOzf "+filename_str+' '+PDBfile);
            while(fp_gz2.good())
            {
                getline(fp_gz2,line);
                if (line.length()<53||line.substr(0,3)=="END") continue;
                parse_pdb_line(line,pep,chain,residue,atom,PDBmap[PDBfile],
                    atomic_detail,allowX);
            }
            fp_gz2.close();
        }

        chain.residues.clear();
        residue.atoms.clear();
        chainIDmap.clear();
        PDBmap.clear();
        return pep;
    }
    

    int use_stdin=(filename_str=="-");
    int use_pstream=0; // input is compressed

    ifstream fp;
    redi::ipstream fp_gz; // if file is compressed
    if (filename_str.length()>=3 && 
        filename_str.substr(filename_str.length()-3,3)==".gz")
    {
        // gzip pdb
        fp_gz.open("zcat "+filename_str);
        use_pstream=1;
    }
    else
    {
        fp.open(filename,ios::in); //ifstream fp(filename,ios::in);
    }

    while(use_stdin?cin.good():(use_pstream?fp_gz.good():fp.good()))
    {
        if (use_stdin)
            getline(cin,line);
        else if (use_pstream)
            getline(fp_gz,line);
        else
            getline(fp,line);

        if (line.length()<53) continue;
        if (line.substr(0,3)=="END") break;
        
        parse_pdb_line(line,pep,chain,residue,atom,chainIDmap,
            atomic_detail,allowX);
    }
    if (!use_stdin)
    {
        if (use_pstream==0)
            fp.close();
        else
            fp_gz.close();
    }
    
    chain.residues.clear();
    residue.atoms.clear();
    chainIDmap.clear();
    return pep;
}

/* filename - full output filename, write to stdout if filename=="-" */
void write_pdb_structure(const char *filename,ModelUnit &pep)
{
    stringstream buf;

    AtomUnit atom;
    ResidueUnit residue;
    ChainUnit chain;
    int a,r,c,i=1;
    for (c=0;c<pep.chains.size();c++)
    {
        chain=pep.chains[c];
        for (r=0;r<chain.residues.size();r++)
        {
            residue=chain.residues[r];
            for (a=0;a<residue.atoms.size();a++)
            {
                atom=residue.atoms[a];
                buf<<setiosflags(ios::left)<<setw(6)
                   <<(residue.het?"HETATM":"ATOM")<<resetiosflags(ios::left)
                   <<setw(5)<<i++<<' '<<atom.name<<' '<<residue.resn<<' '
                   <<chain.chainID<<setw(4)<<residue.resi<<residue.icode
                   <<"   " <<setiosflags(ios::fixed)<<setprecision(3)
                   <<setw(8)<<atom.xyz[0]<<setw(8)<<atom.xyz[1]
                   <<setw(8)<<atom.xyz[2]<<endl;
            }
        }
        buf<<"TER"<<endl;
    }
    
    if (strcmp(filename,"-")==0)
        cout<<buf.str();
    else
    {
        ofstream fp(filename);
        fp<<buf.str();
        fp.close();
    }
}

void write_pdb_structure(const char *filename,ChainUnit &chain)
{
    ModelUnit pep;
    pep.chains.push_back(chain);
    write_pdb_structure(filename,pep);
    pep.chains.clear();
}

/* renumber residue index from "startindex" */
void reindex_pdb(const int startindex,ChainUnit& chain)
{
    for (int r=0;r<chain.residues.size();r++)
        chain.residues[r].resi=r+startindex;
}

void reindex_pdb(const int startindex,ModelUnit& pep)
{
    for (int c=0;c<pep.chains.size();c++) 
        reindex_pdb(startindex,pep.chains[c]);
}

/* convert pdb structure to fasta sequence 
 * convertX - how to deal with non-standard amino acids
 *            0 - only 20 standard amino acids
 *            1 - 20 standard amino acids + MSE
 *            2 - non-standard amino acid with known parent, 
 *                all to legal amino acid in BLOSUM
 *            3 - non-standard amino acid with known parent
 */
inline char aa3to1(const string resn,const int convertX=2)
{
    // 20 standard amino acid + MSE
    if (resn=="ALA") return 'A';
    if (resn=="CYS") return 'C';
    if (resn=="ASP") return 'D';
    if (resn=="GLU") return 'E';
    if (resn=="PHE") return 'F';
    if (resn=="GLY") return 'G';
    if (resn=="HIS") return 'H';
    if (resn=="ILE") return 'I';
    if (resn=="LYS") return 'K';
    if (resn=="LEU") return 'L';
    if (resn=="MET") return 'M';
    if (resn=="ASN") return 'N';
    if (resn=="PRO") return 'P';
    if (resn=="GLN") return 'Q';
    if (resn=="ARG") return 'R';
    if (resn=="SER") return 'S';
    if (resn=="THR") return 'T';
    if (resn=="VAL") return 'V'; 
    if (resn=="TRP") return 'W';
    if (resn=="TYR") return 'Y';

    if (resn=="MSE" && convertX>=1) return 'M';

    if (convertX>=2)
    {
        // non-standard amino acid with known parent
        if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
            resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
            resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
            return 'A';
        if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
            resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
            resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
            resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
            resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
            resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
        if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
            resn=="PCA") return 'E';
        if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
            resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
        if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
            return 'F';
        if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
            resn=="MPQ"||resn=="NMC") return 'G';
        if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
            resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
            return 'H';
        if (resn=="IIL"||resn=="DIL") return 'I';
        if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
            resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
        if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
            resn=="NLN"||resn=="MLE") return 'L';
        if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
        if (resn=="MEN") return 'N';
        if (resn=="DPR"||resn=="HYP") return 'P';
        if (resn=="DGN") return 'Q';
        if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
            resn=="ARM") return 'R';
        if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
            resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
        if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
        if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
        if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
            return 'W';
        if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
            resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
        
        // undeterminted amino acid
        if (resn=="ASX") return 'B'; // or D or N
        if (resn=="GLX") return 'Z'; // or Q or E

        // amino acid with no code in BLOSUM62
        if (convertX>=3)
        {
            if (resn=="SEC") return 'U';
            if (resn=="PYL") return 'O';
        }
        if (resn=="SEC") return 'C';
        if (resn=="PYL") return 'K';
    }
    return 'X';
}

string pdb2fasta(ChainUnit& chain)
{
    chain.sequence="";
    for (int r=0;r<chain.residues.size();r++)
        chain.sequence+=aa3to1(chain.residues[r].resn);
    return chain.sequence;
}

/* ShowSeqLen - whether to show residue number for each chain */
string pdb2fasta(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=0)
{
    stringstream buf;
    string sequence="";
    for (int c=0;c<pep.chains.size();c++)
    {
        sequence=pdb2fasta(pep.chains[c]);
        buf<<'>'<<PDBid<<':'<<pep.chains[c].chainID_full;
        if (ShowSeqLen) buf<<'\t'<<sequence.length();
        buf<<'\n'<<sequence<<'\n';
    }
    sequence.clear();
    return buf.str();
}

/* BackBoneTorsion - backbone torsion angles (Omega, Phi, Psi) */
vector<vector<float> > BackBoneTorsion(ChainUnit& chain)
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
