const char* docstring=""
"qTMclust list PDB/\n"
"    use TM-score to cluster PDB files under folder 'PDB/' and listed by 'list'\n"
"\n"
"input:\n"
"    list - text file listing PDB files for clustering\n"
"    PDB/ - folder for input PDB coordinate files\n"
"\n"
"output:\n"
"    seq.fasta   - fasta format file for amino acid sequence\n"
"    sarst.fasta - fasta format file for SARST code\n"
"    cluster.txt - index file for clustering result\n"
"    ca.xyz      - xyz format file for CA atoms of representative PDB files\n"
"    TM_fast.txt - matrix of TM-score by fast TMalign\n"
"    TM_full.txt - matrix of TM-score by standard TMalign\n"
"\n"
"options:\n"
"    -norm={0,1} protein length with which TM-score is normalized\n"
"        0 - (default) use the larger protein length for normalization\n"
"        1 - use the smaller protein length for normalization\n"
"    -TMmin=0.50     minimum TM-score to consider\n"
"    -TMmax=0.80     maximum TM-score to consider\n"
"    -TMstep=0.10    step size of TM-score cut-offs\n"
"    -MinClustSize=2 minimum cluster size below which clustering will not\n"
"                    be performed\n"
"    -CacheCoor={-1,0,1} whether cache all PDB coordinate in memory\n"
"        0 - (default) cache all coordinate if less than 10,000 chains\n"
"        1 - always cache all coordinates\n"
"       -1 - never cache all coordinates\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "qTMclust.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    int norm=0; // use longer protein
    int MinClustSize=2; // minimum size of cluster
    double TMmin=0.5;  // minimum TM-score cutoff
    double TMmax=0.8;  // maximum TM-score cutoff
    double TMstep=0.1; // step size of TM-score cutoffs
    int CacheCoor=0;   // decide whether cache all PDB based on pdb_entry_num
    vector<string> argv_list;
    for (int arg=1;arg<argc;arg++)
    {
        if (string(argv[arg]).substr(0,6)=="-norm=")
            norm=atoi(string(argv[arg]).substr(6).c_str());
        else if (string(argv[arg]).substr(0,7)=="-TMmin=")
            TMmin=atof(string(argv[arg]).substr(7).c_str());
        else if (string(argv[arg]).substr(0,7)=="-TMmax=")
            TMmax=atof(string(argv[arg]).substr(7).c_str());
        else if (string(argv[arg]).substr(0,8)=="-TMstep=")
            TMstep=atof(string(argv[arg]).substr(8).c_str());
        else if (string(argv[arg]).substr(0,14)=="-MinClustSize=")
            MinClustSize=atof(string(argv[arg]).substr(14).c_str());
        else if (string(argv[arg]).substr(0,6)=="-CacheCoor=")
            CacheCoor=atoi(string(argv[arg]).substr(11).c_str());
        else
            argv_list.push_back(argv[arg]);
    }
    if(argv_list.size()<2)
    {
        cerr<<docstring;
        return 0;
    }

    /* parse filename */
    vector<string> pdb_name_list;
    vector<string> pdb_file_list;
    int pdb_entry_num=parse_pdb_list(argv_list[0],argv_list[1],
        pdb_name_list,pdb_file_list);
    if (CacheCoor==0) CacheCoor=2*(pdb_entry_num<=10000)-1;

    /* parse pdb_file */
    cout<<"convert "<<pdb_entry_num<<" chains into SARST sequence"<<endl;
    vector<pair<int,ChainUnit> > pdb_chain_list;
    vector<pair<int,string> > pdb_name_pair_list;
    vector<pair<int,string> > pdb_file_pair_list;
    ModelUnit tmp_model;
    int i,j;
    int L;
    for (i=0;i<pdb_entry_num;i++)
    {
        cout<<"    parse "<<pdb_file_list[i]<<'\t'<<setiosflags(ios::fixed)
            <<setprecision(2)<<100.*(i+1)/pdb_entry_num<<'%'<<endl;
        // backbone; MSE to MET
        tmp_model=read_pdb_structure(pdb_file_list[i].c_str(),1,1);
        pdb2fasta(tmp_model.chains[0]);
        pdb2sarst(tmp_model.chains[0]);
        
        //free memory for caching coordinates, which are re-parsed by qTMclust
        if (CacheCoor==-1) tmp_model.chains[0].residues.clear();
        else remove_sidechain(tmp_model.chains[0],0); // remove non-CA atoms

        L=tmp_model.chains[0].sarst.length();
        pdb_chain_list.push_back(make_pair(L,tmp_model.chains[0]));
        pdb_name_pair_list.push_back(make_pair(L,pdb_name_list[i]));
        pdb_file_pair_list.push_back(make_pair(L,pdb_file_list[i]));
        tmp_model.chains.clear();
    }

    /* sort pdb list in descending order of chain length */
    stable_sort(pdb_chain_list.begin(),pdb_chain_list.end(),
        cf_pdb_chain_list);
    stable_sort(pdb_name_pair_list.begin(),pdb_name_pair_list.end(),
        cf_int_str_pair_list);
    stable_sort(pdb_file_pair_list.begin(),pdb_file_pair_list.end(),
        cf_int_str_pair_list);
    for (i=0;i<pdb_entry_num;i++)
    {
        pdb_name_list[i]=pdb_name_pair_list[i].second;
        pdb_file_list[i]=pdb_file_pair_list[i].second;
    }
    pdb_name_pair_list.clear();
    pdb_file_pair_list.clear();


    /* write fasta */
    cout<<"write sarst.fasta"<<endl;
    ofstream fp_sarst("sarst.fasta");
    for (i=0;i<pdb_entry_num;i++) fp_sarst<<'>'<<pdb_name_list[i]<<'\t'<<
        pdb_chain_list[i].first<<endl<<pdb_chain_list[i].second.sarst<<endl;
    fp_sarst.close();
    cout<<"write seq.fasta"<<endl;
    ofstream fp_aa("seq.fasta");
    for (i=0;i<pdb_entry_num;i++) fp_aa<<'>'<<pdb_name_list[i]<<'\t'<<
        pdb_chain_list[i].first<<endl<<pdb_chain_list[i].second.sequence<<endl;
    fp_aa.close();

    /* matrix for TM-score */
    cout<<"allocate TM-score table"<<endl;
    // (i,j) store TM-score between i and j, as normalized by i
    // tm_fast_mat is calculated by TMalign -f 8
    // tm_full_mat is calculated by standard TMalign
    // use unsigned char instead of float or double save 4x ~ 8x space
    vector <unsigned char> tmp_array(pdb_entry_num,0);
    vector<vector<unsigned char> >tm_fast_mat(pdb_entry_num,tmp_array);
    map<int,map<int,unsigned char> >tm_full_mat; // sparse matrix
    //vector<vector<unsigned char> >tm_full_mat(pdb_entry_num,tmp_array);
    tmp_array.clear();
    TMclustUnit TMclust;
    initialize_TMclust(TMclust,pdb_entry_num);

    /* perform initial clustering */
    cout<<"heuristic clustering for TM-score "<<TMmin<<endl;
    qTMclust(TMclust, pdb_name_list, pdb_file_list, 
        tm_fast_mat, tm_full_mat, pdb_chain_list, TMmin, 8, norm);
    //fast_clustering(TMclust, pdb_name_list, 
        //tm_fast_mat, tm_full_mat, pdb_chain_list, TMmin, 8, norm);

    /* output initial clusters */
    cout<<"write output for TM-score "<<TMmin<<endl;
    write_TMclust_result("cluster.txt",TMclust,pdb_name_list,TMmin);
    write_matrix("TM_fast.txt",tm_fast_mat);
    write_matrix("TM_full.txt",tm_full_mat,pdb_entry_num);
    write_TMclust_ca_xyz("ca.xyz", TMclust.repr_list, 
        pdb_name_list, pdb_chain_list);

    if (TMmin+TMstep>=TMmax) return 0;

    /* clustering within clusters */
    full_clustering(TMclust, pdb_name_list, pdb_file_list, pdb_chain_list, 
        tm_fast_mat, tm_full_mat, TMmin, TMmax, TMstep, "cluster.txt",
        'a', norm, MinClustSize, CacheCoor);
    write_matrix("TM_full.txt",tm_full_mat,pdb_entry_num);

    /* clean up */
    cout<<"finished clustering of "<<pdb_entry_num<<" chains into "
        <<TMclust.repr_list.size()<<" clusters."<<endl;
    pdb_chain_list.clear();
    TMclust.repr_list.clear();
    TMclust.clust_list.clear();
    return 0;
}
