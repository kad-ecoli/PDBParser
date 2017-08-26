/* get rotamer sequence based on chi-1 angle */
#include <cstring>
#include "SidechainTorsion.hpp"

using namespace std;

char getRotSeq(ResidueUnit& residue,double chi1_angle)
{
    if (chi1_angle<0) chi1_angle+=360;
    if (residue.resn=="ALA")                         return 'A';
    if (residue.resn=="CYS")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'B';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'C';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'D';
    if (residue.resn=="ASP")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'E';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'F';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'G';
    if (residue.resn=="GLU")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'H';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'I';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'J';
    if (residue.resn=="PHE")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'K';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'L';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'M';
    if (residue.resn=="GLY")                         return 'N';
    if (residue.resn=="HIS")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'O';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'P';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'Q';
    if (residue.resn=="ILE")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'R';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'S';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'T';
    if (residue.resn=="LYS")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'U';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'V';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'W';
    if (residue.resn=="LEU")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'X';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'Y';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'Z';
    if (residue.resn=="MET")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'a';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'b';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'c';
    if (residue.resn=="ASN")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'd';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'e';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'f';
    if (residue.resn=="PRO") // perhaps 180 ?
        if      (chi1_angle>=  0 && chi1_angle< 240) return 'g';
        else if (chi1_angle>=240 && chi1_angle< 360) return 'h';
        else if (chi1_angle==360) return 'g';
    if (residue.resn=="GLN")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'i';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'j';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'k';
    if (residue.resn=="ARG")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'l';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'm';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'n';
    if (residue.resn=="SER")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'o';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'p';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'q';
    if (residue.resn=="THR")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'r';
        else if (chi1_angle>=120 && chi1_angle< 240) return 's';
        else if (chi1_angle>=240 && chi1_angle<=360) return 't';
    if (residue.resn=="VAL")
        if      (chi1_angle>=  0 && chi1_angle<120) return 'u';
        else if (chi1_angle>=120 && chi1_angle<240) return 'v';
        else if (chi1_angle>=240 && chi1_angle<360) return 'w';
        else if (chi1_angle==360) return 'v';
    if (residue.resn=="TRP")
        if      (chi1_angle>=  0 && chi1_angle< 120) return 'x';
        else if (chi1_angle>=120 && chi1_angle< 240) return 'y';
        else if (chi1_angle>=240 && chi1_angle<=360) return 'z';
    if (residue.resn=="TYR")
        if      (chi1_angle>=  0 && chi1_angle< 120) return '0';
        else if (chi1_angle>=120 && chi1_angle< 240) return '1';
        else if (chi1_angle>=240 && chi1_angle<=360) return '2';
    return '3';
}

string getRotSeq(ChainUnit& chain)
{
    string RotSeq="";
    bool chi1_only=true;
    vector<vector<double> >angle_mat=SidechainTorsion(chain,chi1_only);
    for (int r=0;r<chain.residues.size();r++)
        if (has_atom_name(chain.residues[r]))
            RotSeq+=getRotSeq(chain.residues[r],angle_mat[r][0]);
    return RotSeq;
}
