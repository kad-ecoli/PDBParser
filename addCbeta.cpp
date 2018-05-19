const char* docstring=
"addCbeta backbone.pdb beta.pdb\n"
"    add C beta atom to residues that lack C beta\n"
"\n"
"Options:\n"
"    -gly_beta={true,false} whether to add C beta on gly\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "PDBParser.hpp"
#include "PeptideBuilder.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    bool gly_beta=true;
    vector<string> argvs;
    string arg;
    for (int a=1;a<argc;a++)
    {
        arg=string(argv[a]);
        if (arg.substr(0,10)=="-gly_beta=")
            gly_beta=(arg.substr(10)=="true");
        else if (arg[0]=='-')
        {
            cerr<<"ERROR! No such option "<<arg<<endl;
            return 0;
        }
        else argvs.push_back(arg);
        arg.clear();
    }

    if(argvs.size()==0)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=2; // read full atom
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argvs[0].c_str(),
        atomic_detail,allowX);

    /* add C beta */
    Geo geo;
    addCbeta(pdb_entry,geo,gly_beta);

    /* make structure */
    write_pdb_structure((argvs.size()>1?argvs[1].c_str():"-"),pdb_entry);

    /* clean up */
    pdb_entry.chains.clear();
    argvs.clear();
    return 0;
}
