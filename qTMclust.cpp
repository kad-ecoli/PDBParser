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

    /* parse pdb_file */
    vector<pair<int,ChainUnit> > pdb_chain_list;
    vector<pair<int,string> > pdb_name_pair_list;
    vector<pair<int,string> > pdb_file_pair_list;
    ChainUnit tmp_chain;
    int i,j;
    for (i=0;i<pdb_entry_num;i++)
    {
        cout<<"    parsing "<<pdb_file_list[i]<<endl;
        // backbone; MSE to MET
        tmp_chain=read_pdb_structure(pdb_file_list[i].c_str(),1,1).chains[0];
        pdb2fasta(tmp_chain);
        pdb2sarst(tmp_chain);
        remove_sidechain(tmp_chain,0); // remove non-CA atoms

        pdb_chain_list.push_back(make_pair(
            tmp_chain.residues.size(),tmp_chain));
        pdb_name_pair_list.push_back(make_pair(
            tmp_chain.residues.size(),pdb_name_list[i]));
        pdb_file_pair_list.push_back(make_pair(
            tmp_chain.residues.size(),pdb_file_list[i]));
    }

    /* sort pdb list in descending order of chain length */
    sort(pdb_chain_list.begin(),pdb_chain_list.end(),
        cf_pdb_chain_list);
    sort(pdb_name_pair_list.begin(),pdb_name_pair_list.end(),
        cf_int_str_pair_list);
    sort(pdb_file_pair_list.begin(),pdb_file_pair_list.end(),
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

    /* convert SARST to int */
    vector<vector<int> >sarst2int_list;
    for (i=0;i<pdb_entry_num;i++)
        sarst2int_list.push_back(sarst2int(pdb_chain_list[i].second.sarst));

    /* matrix for TM-score */
    vector <double> tmp_array(pdb_entry_num,0);
    // (i,j) store TM-score between i and j, as normalized by i
    vector<vector<double> >tm_fast_mat(pdb_entry_num,tmp_array);//TMalign-f8
    vector<vector<double> >tm_full_mat(pdb_entry_num,tmp_array);//TMalign

    /* NWalign using sarst */
    cout<<"SARST alignment derived RMSD superposition"<<endl;
    batch_sarst_rmsd(pdb_name_list,pdb_chain_list, sarst2int_list,
        tm_fast_mat, tm_full_mat, TMmin, norm);
    sarst2int_list.clear();
    write_matrix("TM_fast.txt",tm_fast_mat);
    write_matrix("TM_full.txt",tm_full_mat);

    /* perform initial clustering */
    TMclustUnit TMclust;
    initialize_TMclust(TMclust,pdb_entry_num);
    cout<<"heuristic clustering for TM-score "<<TMmin<<endl;
    fast_clustering(TMclust, pdb_name_list, 
        tm_fast_mat, tm_full_mat, pdb_chain_list, TMmin, 8, norm);

    /* output initial clusters */
    cout<<"write output for TM-score "<<TMmin<<endl;
    write_TMclust_result("cluster.txt",TMclust,pdb_name_list,TMmin);
    write_matrix("TM_fast.txt",tm_fast_mat);
    write_matrix("TM_full.txt",tm_full_mat);
    write_TMclust_ca_xyz("ca.xyz", TMclust.repr_list, 
        pdb_name_list, pdb_chain_list);

    if (TMmin+TMstep>=TMmax) return 0;

    /* clustering within clusters */
    string clust_txt=full_clustering(TMclust, pdb_name_list, pdb_chain_list, 
        tm_fast_mat, tm_full_mat, TMmin, TMmax, TMstep, "cluster.txt",
        "ca.xyz", 'a', norm, MinClustSize);
    write_matrix("TM_full.txt",tm_full_mat);
    return 0;
}
