const char* docstring=
"PeptideBuilder seq.fasta rama.txt out.pdb\n"
"    reconstruct three-dimensional structure for query sequence 'seq.fasta'\n"
"    using backbone torsion angle table 'rama.txt', and output structure\n"
"    to 'out.pdb'\n"
"\n"
"    rama.txt is a three-column table, with each column specifying the\n"
"    omega, phi, psi torsion angles in degrees for one residue.\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "PDBParser.hpp"
#include "PeptideBuilder.hpp"

using namespace std;

string readSingleSequenceFasta(const char* filename)
{   
    ifstream fp(filename, ios::in);
    int use_stdin=(strcmp(filename,"-")==0);
    if (!fp && !use_stdin)
    {
        cerr<<"ERROR! Cannot read file "<<filename<<endl;
        return 0;
    }

    string line,sequence;
    while (use_stdin?cin.good():fp.good())
    {
        use_stdin?getline(cin,line):getline(fp,line);
        if (!line.empty() && line[0]!='>') sequence+=line;
    }
    return sequence;
}

void readRamaTable(const char* filename, vector<vector<float> >&rama_table)
{
    ifstream fp(filename, ios::in);
    int use_stdin=(strcmp(filename,"-")==0);
    if (!fp && !use_stdin)
    {
        cerr<<"ERROR! Cannot read file "<<filename<<endl;
        return;
    }

    string line,sequence;
    int i=0;
    while (use_stdin?cin.good():fp.good())
    {
        if (use_stdin) 
            cin>>rama_table[i][0]>>rama_table[i][1]>>rama_table[i][2];
        else fp>>rama_table[i][0]>>rama_table[i][1]>>rama_table[i][2];
        if (++i >= rama_table.size()) break;
    }
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }

    /* read query sequence */
    string sequence=readSingleSequenceFasta(argv[1]);
    int L=sequence.size();

    /* read backbone torsion table */
    vector<float> tmp_array(3,180);
    vector<vector<float> >rama_table(L,tmp_array);
    readRamaTable(argv[2],rama_table);
    
    /* make structure */
    ChainUnit chain;
    make_chain(chain, sequence, rama_table);    
    write_pdb_structure((argc>3?argv[3]:"-"),chain);

    /* clean up */
    tmp_array.clear();
    rama_table.clear();
    sequence.clear();
    return 0;
}
