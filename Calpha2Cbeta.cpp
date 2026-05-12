const char* docstring=
"Calpha2Cbeta2 ca.pdb beta.pdb\n"
"    add C beta atom to C alpha trace\n"
;

#include <iostream>
#include <vector>
#include <string>

#include "PDBParser.hpp"
#include "GeometryTools.hpp"
#include "StructuralAlphabet.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    vector<string> argvs;
    string arg;
    for (int a=1;a<argc;a++)
    {
        arg=string(argv[a]);
        argvs.push_back(arg);
        arg.clear();
    }

    if(argvs.size()==0)
    {
        cerr<<docstring;
        return 0;
    }

    int atomic_detail=0; // only read CA atom
    int allowX=1;        // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(argvs[0].c_str(),
        atomic_detail,allowX);
    AtomUnit atom;
    atom.name=" CB ";
    atom.xyz.assign(3,0);
    float CACB=1.53;
    vector<float> v1(3,0);
    vector<float> v2(3,0);
    vector<float> v3(3,0);
    int L;

    /* add C beta */
    stringstream buf;
    for (int c=0;c<pdb_entry.chains.size();c++)
    {
        L=pdb_entry.chains[c].residues.size();
        if (L<=2)
        {
            if (L==1 && !(pdb_entry.chains[c].residues[0].resn=="GLY"))
            {
                atom.xyz[0]=pdb_entry.chains[c].residues[0].atoms[0].xyz[0];
                atom.xyz[1]=pdb_entry.chains[c].residues[0].atoms[0].xyz[1];
                atom.xyz[2]=pdb_entry.chains[c].residues[0].atoms[0].xyz[2]
                           +CACB;
                pdb_entry.chains[c].residues[0].atoms.push_back(atom);
            }
            else if (L==2)
            {
                subtract(
                    pdb_entry.chains[c].residues[0].atoms[0].xyz,
                    pdb_entry.chains[c].residues[1].atoms[0].xyz,
                    v1);
                norm(v1,v2);
                v2[0]*=CACB;
                v2[1]*=CACB;
                v2[2]*=CACB;
                atom.xyz[0]=pdb_entry.chains[c].residues[0].atoms[0].xyz[0]+v2[0];
                atom.xyz[1]=pdb_entry.chains[c].residues[0].atoms[0].xyz[1]+v2[1];
                atom.xyz[2]=pdb_entry.chains[c].residues[0].atoms[0].xyz[2]+v2[2];
                if (!(pdb_entry.chains[c].residues[0].resn=="GLY"))
                    pdb_entry.chains[c].residues[0].atoms.push_back(atom);
                atom.xyz[0]=pdb_entry.chains[c].residues[1].atoms[0].xyz[0]-v2[0];
                atom.xyz[1]=pdb_entry.chains[c].residues[1].atoms[0].xyz[1]-v2[1];
                atom.xyz[2]=pdb_entry.chains[c].residues[1].atoms[0].xyz[2]-v2[2];
                if (!(pdb_entry.chains[c].residues[1].resn=="GLY"))
                    pdb_entry.chains[c].residues[1].atoms.push_back(atom);
            }
            continue;
        }

        pdb2ss(pdb_entry.chains[c]);
        for (int r=0;r<L;r++)
        {
            if (pdb_entry.chains[c].residues[r].resn=="GLY") continue;
            if (r==0) subtract(
                pdb_entry.chains[c].residues[2].atoms[0].xyz,
                pdb_entry.chains[c].residues[1].atoms[0].xyz,
                v1);
            else subtract(
                pdb_entry.chains[c].residues[r].atoms[0].xyz,
                pdb_entry.chains[c].residues[r-1].atoms[0].xyz,
                v1);
            norm_warnless(v1);
            if (r+1==L) subtract(
                pdb_entry.chains[c].residues[r-2].atoms[0].xyz,
                pdb_entry.chains[c].residues[r-1].atoms[0].xyz,
                v2);
            else subtract(
                pdb_entry.chains[c].residues[r].atoms[0].xyz,
                pdb_entry.chains[c].residues[r+1].atoms[0].xyz,
                v2);
            norm_warnless(v2);
            vectorsum(v1,v2,v3);
            norm_warnless(v3);
            atom.xyz[0]=pdb_entry.chains[c].residues[r].atoms[0].xyz[0]+v3[0]*CACB;
            atom.xyz[1]=pdb_entry.chains[c].residues[r].atoms[0].xyz[1]+v3[1]*CACB;
            atom.xyz[2]=pdb_entry.chains[c].residues[r].atoms[0].xyz[2]+v3[2]*CACB;

            /* rotate by 20.9 ~ 23.2 ~ 51.2 */
            if (r==0) subtract(
                pdb_entry.chains[c].residues[2].atoms[0].xyz,
                pdb_entry.chains[c].residues[0].atoms[0].xyz,
                v2);
            else if (r+1==L) subtract(
                pdb_entry.chains[c].residues[L-1].atoms[0].xyz,
                pdb_entry.chains[c].residues[L-3].atoms[0].xyz,
                v2);
            else subtract(
                pdb_entry.chains[c].residues[r+1].atoms[0].xyz,
                pdb_entry.chains[c].residues[r-1].atoms[0].xyz,
                v2);
            vectorsum(pdb_entry.chains[c].residues[r].atoms[0].xyz,v2,v3);
            
            //#mean         short    medm    long     all
            //CB           0.9881  0.9797  0.9735  0.9881
            //CA8          0.7395  0.7007  0.6717  0.7068
            //CA9          0.6216  0.6528  0.6412  0.6463
            //pseudoCB(0)  0.8531  0.8381  0.8167  0.8411
            //pseudoCB(20) 0.9007  0.8922  0.8796  0.8986
            //pseudoCB(25) 0.9098  0.9021  0.8924  0.9103
            //pseudoCB(30) 0.9153  0.9075  0.9004  0.9175
            //pseudoCB(35) 0.9155  0.9084  0.9040  0.9204
            //pseudoCB(40) 0.9007  0.8922  0.8796  0.8986
            //pseudoCB(50) 0.8938  0.8862  0.8853  0.9026
            //pseudoCB(20) 0.9328  0.9258  0.9204  0.9380 * // H=51.2, E=23.2, C=35
            //pseudoCB(25) 0.9323  0.9251  0.9198  0.9375   // H=51.2, E=20.9, C=35
            //pseudoCB(30) 0.9327  0.9256  0.9202  0.9378   // H=51.2, E=22.0, C=35
            //pseudoCB(20) 0.9238  0.9174  0.9107  0.9291   // H=51.2, E=23.2, C=20
            //pseudoCB(25) 0.9291  0.9223  0.9163  0.9342   // H=51.2, E=23.2, C=25
            //pseudoCB(30) 0.9323  0.9250  0.9195  0.9374   // H=51.2, E=23.2, C=30
            //pseudoCB(40) 0.9311  0.9244  0.9188  0.9366   // H=51.2, E=23.2, C=40
            //pseudoCB(45) 0.9271  0.9208  0.9148  0.9331   // H=51.2, E=23.2, C=45
            float angle=35;
            switch (pdb_entry.chains[c].ss[r])
            {
                case 'H': angle=51.2; break; // helix
                case 'E': angle=23.2; break; // parallel/antiparallel 20.9/23.2
                case 'S': angle=23.2; break; // parallel/antiparallel 20.9/23.2
                default:  angle=35;
            }

            CoordinateRotation(atom.xyz,
                v3,pdb_entry.chains[c].residues[r].atoms[0].xyz,
                angle,v1);
            atom.xyz[0]=v1[0];
            atom.xyz[1]=v1[1];
            atom.xyz[2]=v1[2];

            pdb_entry.chains[c].residues[r].atoms.push_back(atom);
        }
        //vector<bool> tmp_vec(L,0);
        //vector<vector<bool> > ct_mat(L,tmp_vec);
        for (int r1=0;r1<L;r1++)
        {
            for (int r2=r1+6;r2<L;r2++)
            {
                if (Points2Distance(
                    pdb_entry.chains[0].residues[r1].atoms.back().xyz,
                    pdb_entry.chains[0].residues[r2].atoms.back().xyz)>8)
                    continue;
                //ct_mat[r1][r2]=ct_mat[r2][r1]=1;
                //if (abs(r2-r1)<6) continue;
                buf<<r1+1<<' '<<r2+1<<" 0 8 1"<<endl;
            }
        }
        //vector<bool>().swap(tmp_vec);
        //vector<vector<bool> > ().swap(ct_mat);
    }

    /* make structure */
    write_pdb_structure((argvs.size()>1?argvs[1].c_str():"-"),pdb_entry);

    if (argvs.size()>2)
    {
        ofstream fp(argvs[2]);
        fp<<buf.str();
        fp.close();
    }
    buf.str(string());

    /* clean up */
    pdb_entry.chains.clear();
    argvs.clear();
    return 0;
}
