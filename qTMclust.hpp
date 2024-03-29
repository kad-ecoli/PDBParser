#ifndef QTMCLUST_HPP
#define QTMCLUST_HPP 1

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
#include "FilePathParser.hpp"
#include "PDBParser.hpp"
#include "MathTools.hpp"
#include "NWalign.hpp"
#include "pdb2rmsd.hpp"
#include "Superpose.hpp"
#include "TMalign.hpp"

using namespace std;

struct TMclustUnit // struc for storing clustering result
{
    vector<size_t> unclust_list; // index of unclustered entry
    vector<size_t> repr_list;    // index of cluster representative
    vector<vector<size_t> > clust_list; // each row is one cluster
                              // each column is one member in a cluster
    map<size_t,size_t> clust_map;   // key is entry, value is cluster index
};

int parse_pdb_list(const string pdb_list_file,const string pdb_folder,
    vector<string>&pdb_name_list,vector<string>&pdb_file_list)
{
    string pdb_folder_fullname=string(pdb_folder);
    if (pdb_folder_fullname.length()>0 &&
        pdb_folder_fullname[pdb_folder_fullname.length()-1]!='/')
        pdb_folder_fullname+='/';

    ifstream fp(pdb_list_file.c_str(),ios::in);
    string line;
    int nonexist_file_count=0;

    vector <string> prefix_list;
    prefix_list.push_back("");
    prefix_list.push_back("pdb");
    
    vector <string> suffix_list;
    suffix_list.push_back("");
    suffix_list.push_back(".pdb");
    suffix_list.push_back(".ent");
    suffix_list.push_back(".pdb.gz");
    suffix_list.push_back(".ent.gz");
    suffix_list.push_back("-pdb-bundle.tar.gz");
    suffix_list.push_back(".tar.gz");

    int p,s,i;
    bool pdb_file_found=false;
    bool found_dup=false;
    string pdb_file_fullname;
    while(fp.good())
    {
        getline(fp,line);
        if (line.length()==0 || line[0]=='#') continue;

        bool found_dup=false;
        for (i=0;i<pdb_name_list.size();i++)
        {
            if (pdb_name_list[i]==line)
            {
                found_dup=true;
                break;
            }
        }
        if (found_dup)
        {
            cout<<"Warning! skip duplicated entry "<<line<<endl;
            continue;
        }

        bool pdb_file_found=false;
        for (p=0;p<prefix_list.size()&&pdb_file_found==false;p++)
        {
            for (s=0;s<suffix_list.size()&&pdb_file_found==false;s++)
            {
                pdb_file_fullname=pdb_folder_fullname+
                    prefix_list[p]+line+suffix_list[s];
                if (isfile(pdb_file_fullname.c_str()))
                {
                    pdb_file_found=true;
                    break;
                }
            }
        }

        if (pdb_file_found)
        {
            pdb_name_list.push_back(line);
            pdb_file_list.push_back(pdb_file_fullname);
        }
        else cerr<<"ERROR! Cannot locate PDB file for "<<line<<endl;
    }
    fp.close();
    return pdb_name_list.size();
}

/* function for sorting vector<pair<size_t,ChainUnit> > */ 
bool cf_pdb_chain_list(pair<size_t,ChainUnit> chain1,pair<size_t,ChainUnit> chain2)
{
    return chain1.first>chain2.first;
}

bool cf_int_str_pair_list(pair<size_t,string> chain1,pair<size_t,string> chain2)
{
    return chain1.first>chain2.first;
}

/* NWalign+Kabsch rough TM-score calculation */
void rough_tmscore_by_kabsch(const string& seq1, const string& seq2, 
    const vector<int>& seq2int1, const vector<int>& seq2int2,
    string & aln1,string & aln2, ChainUnit& chain1, ChainUnit& chain2,
    vector<vector<float> >&xyz_list1, vector<vector<float> >&xyz_list2, 
    vector<vector<float> >&RotMatix, vector<float>&TranVect,
    float& tmscore1, float& tmscore2, const int L1, const int L2,
    const int seq_type=0, const int glocal=0)
{
    NWalign(seq1, seq2, seq2int1, seq2int2, aln1, aln2, seq_type, glocal);
    aln2coor(aln1, aln2, chain1, chain2, xyz_list1, xyz_list2,0);
    int aln_len=xyz_list1.size();
    if (aln_len==0) tmscore1=tmscore2=0;
    else
    {
        /* kabsch */
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

        /* change coordinate */
        vector<float> tmp_array(3,0);
        vector<vector<float> > super_xyz_list1(aln_len,tmp_array);
        for (int r=0; r<aln_len; r++)
            ChangeCoor(xyz_list1[r], RotMatix, TranVect, super_xyz_list1[r]);

        /* TM-score */
        tmscore1=calTMscore(super_xyz_list1,xyz_list2,L1);
        tmscore2=calTMscore(super_xyz_list1,xyz_list2,L2);

        /* clean-up */
        super_xyz_list1.clear();
    }
}

void initialize_TMclust(TMclustUnit &TMclust,const size_t pdb_entry_num)
{
    TMclust.repr_list.clear();
    TMclust.clust_list.clear();
    TMclust.unclust_list.clear();

    vector<size_t> tmp_array(1,0);
    TMclust.repr_list.push_back(0); // first chain is first cluster
    TMclust.clust_list.push_back(tmp_array);
    for (size_t i=pdb_entry_num-1;i>0;i--) // vector cannot pop_front. so reverse
        TMclust.unclust_list.push_back(i); // the order of member list
}

/* initialize TMclustUnit with members list */
void initialize_TMclust(TMclustUnit &TMclust, const vector<size_t> &member_list)
{
    TMclust.repr_list.clear();
    TMclust.clust_list.clear();
    TMclust.unclust_list.clear();
    if (member_list.size()==0) return;

    vector<size_t> tmp_array(1,member_list[0]);
    TMclust.repr_list.push_back(member_list[0]); // first chain is first cluster
    TMclust.clust_list.push_back(tmp_array);
    for (size_t i=member_list.size()-1;i>0;i--)
        TMclust.unclust_list.push_back(member_list[i]);
}

/* add new_member to cluster represented by clust_repr */
size_t add_chain_to_clust(TMclustUnit &TMclust, const size_t new_member,
    const size_t clust_repr)
{
    size_t k; // cluster index
    for (k=0;k<TMclust.repr_list.size();k++)
    {
        if (TMclust.repr_list[k]==clust_repr)
        {
            TMclust.clust_list[k].push_back(new_member);
            TMclust.clust_map[new_member]=k;
            break;
        }
    }
    if (k==TMclust.repr_list.size())
    {
        cout<<"FATAL ERROR! No such representative "<<clust_repr<<endl;
        exit(0);
    }
    return k; // index of clust_repr in repr_list
}

/* assign new_repr as new cluster */
size_t assign_chain_as_new_clust(TMclustUnit &TMclust, const size_t new_repr)
{
    TMclust.repr_list.push_back(new_repr);
    vector<size_t> tmp_array(1,new_repr);
    TMclust.clust_list.push_back(tmp_array);
    TMclust.clust_map[new_repr]=TMclust.repr_list.size()-1;
    return TMclust.clust_map[new_repr];
}

int qTMclust(TMclustUnit &TMclust, const vector<string>&pdb_name_list,
    const vector<string>&pdb_file_list,
    //vector<vector<unsigned char> >&tm_fast_mat,
    map<long long,unsigned char>&tm_fast_mat,
    map<long long,unsigned char>&tm_full_mat,
    vector<pair<size_t,ChainUnit> >&pdb_chain_list,
    const float tmscore_cutoff=0.5, const int f=8, const int norm=0,
    const int CacheCoor=1, const int heuristic=0)
{
    long long pdb_entry_num=pdb_chain_list.size();

    size_t i,j,k,l;
    int add_to_clust=-1;
    int max_clust_size=0;
    map<int,unsigned char> tmp_map;

    string aln_i,aln_j;
    float tmscore_i,tmscore_j;  // tmscore from sarst based superposition
    float tmscore_aa_i,tmscore_aa_j; // tmscore from aa based superposition
    float tmscore_ss_i,tmscore_ss_j; // tmscore from ss based superposition
    float tmscore; // the tmscore to consider when clustering
    int L_i,L_j;
    int aln_len;

    vector<vector<float> > xyz_list_i,xyz_list_j;
    vector<vector<float> > RotMatix;  // U
    vector<float> TranVect;  // t
    vector<pair<float,size_t> >aln_order_pair;

    /**** convert SARST to int ****/
    vector<vector<int> >sarst2int_list;
    vector<vector<int> >seq2int_list;
    vector<vector<int> >ss2int_list;
    for (i=0;i<pdb_entry_num;i++)
    {
        sarst2int_list.push_back(sarst2int(pdb_chain_list[i].second.sarst));
        seq2int_list.push_back(aa2int(pdb_chain_list[i].second.sequence));
        ss2int_list.push_back(ss2int(pdb_chain_list[i].second.ss));
    }

    /* read first chain */
    cout<<"    assigning "<<pdb_name_list[0]<<" as first cluster\t"
        <<setprecision(2)<<100./pdb_entry_num<<'%'<<endl;
    ModelUnit tmp_model=read_pdb_structure(pdb_file_list[0].c_str(),0,1);
    if (pdb_chain_list[0].second.residues.size()==0)
        pdb_chain_list[0].second.residues=tmp_model.chains[0].residues;
    tmp_model.chains.clear();

    /**** perform superposition ****/
    while(TMclust.unclust_list.size())
    {
        i=TMclust.unclust_list.back(); // try to add protein i to clust_list
        L_i=pdb_chain_list[i].first;
        cout<<"    assigning "<<pdb_name_list[i];
        TMclust.unclust_list.pop_back();
        add_to_clust=-1; // cluster representative to whom cluster i belongs
    
        /* parse current chain on the fly */
        if (pdb_chain_list[i].second.residues.size()==0)
        {
            tmp_model=read_pdb_structure(pdb_file_list[i].c_str(),0,1);
            pdb_chain_list[i].second.residues=tmp_model.chains[0].residues;
            tmp_model.chains.clear();
        }

        aln_order_pair.clear();
        for (k=0;k<TMclust.repr_list.size();k++)
        {
            j=TMclust.repr_list[k]; // representative for cluster k
            L_j=pdb_chain_list[j].first;
            if (norm==0 &&
               (L_i<tmscore_cutoff*L_j || L_j<tmscore_cutoff*L_i))
                continue; // skip unnecessary superposition

            /*** superpose by sarst ***/
            rough_tmscore_by_kabsch(pdb_chain_list[i].second.sarst,
                pdb_chain_list[j].second.sarst,
                sarst2int_list[i], sarst2int_list[j], aln_i, aln_j,
                pdb_chain_list[i].second, pdb_chain_list[j].second,
                xyz_list_i, xyz_list_j, RotMatix, TranVect,
                tmscore_i, tmscore_j, L_i, L_j, 2,0); //2 -sarst; 0 -global

            /* clean up */
            RotMatix.clear();
            TranVect.clear();
            aln_i.clear();
            aln_j.clear();
            xyz_list_i.clear();
            xyz_list_j.clear();

            /*** superpose by amino acid sequence ***
             * CA-only proteins are often not alignable by SARST */
            rough_tmscore_by_kabsch(pdb_chain_list[i].second.sequence,
                pdb_chain_list[j].second.sequence,
                seq2int_list[i], seq2int_list[j], aln_i, aln_j,
                pdb_chain_list[i].second, pdb_chain_list[j].second,
                xyz_list_i, xyz_list_j, RotMatix, TranVect,
                tmscore_aa_i, tmscore_aa_j, L_i, L_j, 0,0); //0 -aa; 0 -global
            tmscore_i=MAX(tmscore_i,tmscore_aa_i);
            tmscore_j=MAX(tmscore_j,tmscore_aa_j);
            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));

            /* clean up */
            RotMatix.clear();
            TranVect.clear();
            aln_i.clear();
            aln_j.clear();
            xyz_list_i.clear();
            xyz_list_j.clear();

            /*** superpose by secondary structure ***/
            rough_tmscore_by_kabsch(pdb_chain_list[i].second.ss,
                pdb_chain_list[j].second.ss,
                ss2int_list[i], ss2int_list[j], aln_i, aln_j,
                pdb_chain_list[i].second, pdb_chain_list[j].second,
                xyz_list_i, xyz_list_j, RotMatix, TranVect,
                tmscore_ss_i, tmscore_ss_j, L_i, L_j, 4,0); //4 -ss; 0 -global
            tmscore_i=MAX(tmscore_i,tmscore_ss_i);
            tmscore_j=MAX(tmscore_j,tmscore_ss_j);
            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));

            /* clean up */
            RotMatix.clear();
            TranVect.clear();
            aln_i.clear();
            aln_j.clear();
            xyz_list_i.clear();
            xyz_list_j.clear();

            /* the following lines pre-terminate sarst+rmsd alignment if the
             * first significant hit is found. this is not only unhelpful
             * for shortening computational time, but also resulting in
             * assignment to cluster whose TM-score is not the best
            if (tmscore>=tmscore_cutoff*.6)
            {
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, tmscore_i,tmscore_j,0);
                tm_full_mat[i*pdb_entry_num+j]=
                    =tm_fast_mat[i][j]=(int)(255*tmscore_i+.5);
                tm_full_mat[j*pdb_entry_num+i]=
                    =tm_fast_mat[j][i]=(int)(255*tmscore_j+.5);
                tmscore=(norm==0?
                    MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
                if (tmscore>=tmscore_cutoff)
                {
                    add_to_clust=j;
                    break;
                }
            }
            */
            aln_order_pair.push_back(make_pair(tmscore,j));
        }

        /* if we do not pre-terminate sarst+rmsd, add_to_clust==-1
        if (add_to_clust>=0)
        {
            cout<<" to "<<pdb_name_list[j]<<endl;
            add_chain_to_clust(TMclust,i,j);
            aln_order_pair.clear();
            continue;
        }
        */

        if (aln_order_pair.size())
        {
            stable_sort(aln_order_pair.begin(), aln_order_pair.end());
            reverse(aln_order_pair.begin(), aln_order_pair.end());
        }

        size_t max_TMalign_num=MIN(aln_order_pair.size(), size_t(heuristic/L_i));
        for (l=0;l<aln_order_pair.size();l++)
        {
            if (heuristic>0 && l>max_TMalign_num && 
                aln_order_pair[l].first<0.4*tmscore_cutoff) break;
            k=TMclust.clust_map[aln_order_pair[l].second]; // cluster index
            j=TMclust.repr_list[k]; // representative for cluster k

            L_j=pdb_chain_list[j].first;
            if (norm==0 && 
               (L_i<tmscore_cutoff*L_j || L_j<tmscore_cutoff*L_i)) continue;
            
            if (tm_full_mat.count(i*pdb_entry_num+j))
            {
                tmscore_i=tm_full_mat[i*pdb_entry_num+j]/255.;
                tmscore_j=tm_full_mat[j*pdb_entry_num+i]/255.;
            }
            else if (tm_fast_mat.count(i*pdb_entry_num+j))
            {
                tmscore_i=tm_fast_mat[i*pdb_entry_num+j]/255.;
                tmscore_j=tm_fast_mat[j*pdb_entry_num+i]/255.;
            }
            else
            {
                tmscore_i=tmscore_j=0;
                string aln_j; // for unknown reason, aln_j must be re-
                              // allocated 
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, tmscore_i, tmscore_j, f, false);
                aln_j.clear();
                tm_fast_mat[i*pdb_entry_num+j]=(int)(255*tmscore_i+.5);
                tm_fast_mat[j*pdb_entry_num+i]=(int)(255*tmscore_j+.5);
            }
            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
           
            if (f>0 && tmscore>=tmscore_cutoff*.9 && 
                tm_full_mat.count(i*pdb_entry_num+j)==0)
            {
                string aln_j;
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, tmscore_i, tmscore_j, 0, false);
                aln_j.clear();
                tm_full_mat[i*pdb_entry_num+j]=(int)(255*tmscore_i+.5);
                tm_full_mat[j*pdb_entry_num+i]=(int)(255*tmscore_j+.5);
            }
            else if (f==0 && 
                tm_full_mat.count(i*pdb_entry_num+j)==0)
            {
                tm_full_mat[i*pdb_entry_num+j]=(int)(255*tmscore_i+.5);
                tm_full_mat[j*pdb_entry_num+i]=(int)(255*tmscore_j+.5);
            }

            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
            if (tmscore>=tmscore_cutoff)
            {
                add_to_clust=j;
                break;
            }
        }

        if (add_to_clust>=0)
        {
            cout<<" to "<<pdb_name_list[j]<<"\t\t";
            add_chain_to_clust(TMclust,i,j);
            if (CacheCoor==-1) pdb_chain_list[i].second.residues.clear();
        }
        else 
        {
            cout<<" as new cluster\t";
            assign_chain_as_new_clust(TMclust,i);
        }
        cout<<setiosflags(ios::fixed)<<setprecision(2)
            <<100.*(i+1)/pdb_entry_num<<"%\t"
            <<l+(add_to_clust>=0)<<'/'<<aln_order_pair.size()<<'\t'
            <<(aln_order_pair.size()?aln_order_pair[l-(l==
               aln_order_pair.size())].first:0)<<'\t'<<L_i<<endl;
        aln_order_pair.clear();
    }

    /* clean-up */
    seq2int_list.clear();
    sarst2int_list.clear();
    ss2int_list.clear();

    /* return max cluster size */
    for (i=0;i<TMclust.clust_list.size();i++) 
        max_clust_size=MAX(max_clust_size,TMclust.clust_list[i].size());
    return max_clust_size;
}

int fast_clustering(TMclustUnit &TMclust,const vector<string>&pdb_name_list,
    //vector<vector<unsigned char> >&tm_fast_mat,
    map<long long,unsigned char>&tm_fast_mat,
    map<long long,unsigned char>&tm_full_mat,
    vector<pair<size_t,ChainUnit> >&pdb_chain_list,
    const float tmscore_cutoff=0.5, const int norm=0)
{
    long long pdb_entry_num=pdb_chain_list.size();
    size_t i,j,k,l;
    int add_to_clust=-1;
    int max_clust_size=0;
    map<int,unsigned char> tmp_map;

    float tmscore,tmscore_i,tmscore_j;
    int L_i,L_j;

    vector<pair<float,size_t> >aln_order_pair;
    while(TMclust.unclust_list.size())
    {
        i=TMclust.unclust_list.back(); // try to add protein i to clust_list
        L_i=pdb_chain_list[i].first;
        //cout<<"    assigning "<<pdb_name_list[i];
        TMclust.unclust_list.pop_back();
        add_to_clust=-1; // cluster representative to whom cluster i belongs

        /* establish alignment order by tm_fast_mat */
        aln_order_pair.clear();
        for (k=0;k<TMclust.repr_list.size();k++)
        {
            j=TMclust.repr_list[k]; // representative for cluster k
            if (tm_fast_mat.count(i*pdb_entry_num+j))
            {
                tmscore_i=tm_fast_mat[i*pdb_entry_num+j]/255.;
                tmscore_j=tm_fast_mat[j*pdb_entry_num+i]/255.;
            }
            else tmscore_i=tmscore_j=0;
            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
            aln_order_pair.push_back(make_pair(tmscore,j));
        }
        if (aln_order_pair.size())
        {
            stable_sort(aln_order_pair.begin(), aln_order_pair.end());
            reverse(aln_order_pair.begin(), aln_order_pair.end());
        }

        /* add to clust k if full tmscore(i,k)>=tmscore_cutoff */
        for (l=0;l<aln_order_pair.size();l++)
        {
            k=TMclust.clust_map[aln_order_pair[l].second]; // cluster index
            j=TMclust.repr_list[k]; // representative for cluster k

            if (norm==0)
            { 
                L_j=pdb_chain_list[j].first;
                if (L_i<tmscore_cutoff*L_j || L_j<tmscore_cutoff*L_i)
                    continue;
            }

            if (tm_full_mat.count(i*pdb_entry_num+j))
            {
                tmscore_i=tm_full_mat[i*pdb_entry_num+j]/255.;
                tmscore_j=tm_full_mat[j*pdb_entry_num+i]/255.;
            }
            else
            {
                if (tm_fast_mat.count(i*pdb_entry_num+j))
                {
                    tmscore_i=tm_fast_mat[i*pdb_entry_num+j]/255.;
                    tmscore_j=tm_fast_mat[j*pdb_entry_num+i]/255.;
                }
                else tmscore_i=tmscore_j=0;
                tmscore=(norm==0?
                    MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
                if (tmscore>=tmscore_cutoff)
                {
                    add_to_clust=j;
                    break;
                }
                string aln_i,aln_j;
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, tmscore_i,tmscore_j, 0, false);
                aln_i.clear();
                aln_j.clear();
                tm_full_mat[i*pdb_entry_num+j]
                    =tm_fast_mat[i*pdb_entry_num+j]
                    =(int)(255*tmscore_i+.5);
                tm_full_mat[j*pdb_entry_num+i]
                    =tm_fast_mat[j*pdb_entry_num+i]
                    =(int)(255*tmscore_j+.5);
            }
            tmscore=(norm==0?
                MIN(tmscore_i,tmscore_j):MAX(tmscore_i,tmscore_j));
            if (tmscore>=tmscore_cutoff)
            {
                add_to_clust=j;
                break;
            }
        }

        if (add_to_clust>=0)
        {
            add_chain_to_clust(TMclust,i,j);
        }
        else
        {
            TMclust.repr_list.push_back(i);
            vector<size_t> tmp_array(1,i);
            TMclust.clust_list.push_back(tmp_array);
        }
    }
    for (i=0;i<TMclust.clust_list.size();i++) 
        max_clust_size=MAX(max_clust_size,TMclust.clust_list[i].size());
    return max_clust_size;
}

void write_TMclust_result(const string filename,
    const TMclustUnit &TMclust, const vector<string> pdb_name_list,
    float tmscore_cutoff=0.5,char openmode='w')
{
    stringstream buf;
    size_t i,j;

    size_t N=0; // total number of entries
    for (i=0;i<TMclust.clust_list.size();i++) N+=TMclust.clust_list[i].size();
    buf<<"TM_cut="<<setprecision(4)<<tmscore_cutoff
       <<"\tN="<<N<<"\tN_repr="<<TMclust.repr_list.size()<<endl;

    for (i=0;i<TMclust.repr_list.size();i++)
    {
        buf<<pdb_name_list[TMclust.repr_list[i]];
        for (j=1;j<TMclust.clust_list[i].size();j++)
            buf<<'\t'<<pdb_name_list[TMclust.clust_list[i][j]];
        buf<<endl;
    }
    buf<<"$$$$"<<endl;

    ofstream fp;
    if (openmode=='w') fp.open(filename.c_str(),ofstream::trunc);
    else if (openmode=='a') fp.open(filename.c_str(),ofstream::app);
    fp<<buf.str();
    fp.close();
}

void write_matrix(const string filename,
    map<long long,unsigned char>&tm_mat,const vector<string> pdb_name_list)
{
    stringstream buf;
    size_t i,j;
    long long pdb_entry_num=pdb_name_list.size();
    for (i=0;i<pdb_entry_num;i++)
    {
        for (j=1;j<pdb_entry_num;j++)
        {
            if (tm_mat.count(i*pdb_entry_num+j))
            {
                buf<<pdb_name_list[i]<<'\t'<<pdb_name_list[j]<<'\t'
                   <<setprecision(4)<<tm_mat[i*pdb_entry_num+j]/255.<<'\t'
                   <<setprecision(4)<<tm_mat[j*pdb_entry_num+i]/255.<<endl;
            }
        }
    }
    ofstream fp(filename.c_str());
    fp<<buf.str();
    fp.close();
}


void write_matrix(const string filename,
    map<long long,unsigned char>&tm_mat,const long long pdb_entry_num)
{
    stringstream buf;
    size_t i,j;
    for (i=0;i<pdb_entry_num;i++)
    {
        if (tm_mat.count(i*pdb_entry_num)==0) buf<<0;
        else buf<<setprecision(4)<<tm_mat[i*pdb_entry_num]/255.;

        for (j=1;j<pdb_entry_num;j++)
        {
            if (tm_mat.count(i*pdb_entry_num+j)==0) buf<<"\t0";
            else buf<<'\t'<<setprecision(4)<<tm_mat[i*pdb_entry_num+j]/255.;
        }
        buf<<endl;
    }
    ofstream fp(filename.c_str());
    fp<<buf.str();
    fp.close();
}

void write_matrix(const string filename, 
    const vector<vector<unsigned char> >&tm_mat)
{
    stringstream buf;
    size_t i,j;
    for (i=0;i<tm_mat.size();i++)
    {
        buf<<setprecision(4)<<tm_mat[i][0]/255.;
        for (j=1;j<tm_mat[0].size();j++)
        {
            buf<<'\t'<<setprecision(4)<<tm_mat[i][j]/255.;
        }
        buf<<endl;
    }
    ofstream fp(filename.c_str());
    fp<<buf.str();
    fp.close();
}

void write_matrix(const string filename, const vector<vector<float> >&tm_mat)
{
    stringstream buf;
    size_t i,j;
    for (i=0;i<tm_mat.size();i++)
    {
        buf<<setprecision(4)<<tm_mat[i][0];
        for (j=1;j<tm_mat[0].size();j++)
        {
            buf<<'\t'<<setprecision(4)<<tm_mat[i][j];
        }
        buf<<endl;
    }
    ofstream fp(filename.c_str());
    fp<<buf.str();
    fp.close();
}

void write_TMclust_ca_xyz(const string filename, 
    const vector<size_t>&repr_list, const vector<string>& pdb_name_list,
    vector<pair<size_t,ChainUnit> >& pdb_chain_list,char openmode='w')
{
    stringstream buf;
    size_t repr,i,r,res_num;
    for (repr=0;repr<repr_list.size();repr++)
    {
        i=repr_list[repr];
        res_num=0;
        for (r=0;r<pdb_chain_list[i].second.residues.size();r++)
            res_num+=(pdb_chain_list[i].second.residues[r].atoms.size()==1);
        buf<<res_num<<endl<<pdb_name_list[i]<<endl;
        for (r=0;r<pdb_chain_list[i].second.residues.size();r++)
        {
            if (pdb_chain_list[i].second.residues[r].atoms.size()!=1)
                continue;
            buf<<aa3to1(pdb_chain_list[i].second.residues[r].resn)
               <<setiosflags(ios::fixed)<<setprecision(3)
               <<setw(9)<<pdb_chain_list[i].second.residues[r].atoms[0].xyz[0]
               <<setw(9)<<pdb_chain_list[i].second.residues[r].atoms[0].xyz[1]
               <<setw(9)<<pdb_chain_list[i].second.residues[r].atoms[0].xyz[2]
               <<endl;
        }
    }
    ofstream fp;
    if (openmode=='w') fp.open(filename.c_str(),ofstream::trunc);
    else if (openmode=='a') fp.open(filename.c_str(),ofstream::app);
    fp<<buf.str();
    fp.close();
}

int full_clustering(TMclustUnit &TMclust,
    const vector<string>&pdb_name_list, const vector<string>&pdb_file_list, 
    vector<pair<size_t,ChainUnit> >&pdb_chain_list,
    //vector<vector<unsigned char> >&tm_fast_mat,
    map<long long,unsigned char>&tm_fast_mat,
    map<long long,unsigned char>&tm_full_mat,
    float TMmin=0.5, const float TMmax=0.8, const float TMstep=0.1,
    const string cluster_filename="cluster.txt", const char openmode='a',
    const int norm=0, const int MinClustSize=2, const int CacheCoor=1)
{
    TMmin+=TMstep;
    if (TMmin>TMmax) return 0;
    if (TMclust.repr_list.size()==0) return 0;

    cout<<"heuristic clustering for TM-score "<<TMmin<<endl;
    size_t i,j,l;
    int max_clust_size=0;
    ModelUnit tmp_model;
    for (i=0;i<TMclust.repr_list.size();i++)
    {
        if (TMclust.clust_list[i].size()<=MinClustSize) continue;
        cout<<"    subclustering "<<pdb_name_list[TMclust.repr_list[i]]<<endl;

        /* clustering at TMmin */
        TMclustUnit TMsubclust;
        initialize_TMclust(TMsubclust,TMclust.clust_list[i]);
        for (l=0;l<TMclust.clust_list[i].size();l++)
        {
            j=TMclust.clust_list[i][l];
            if (pdb_chain_list[j].second.residues.size()==0)
            {
                tmp_model=read_pdb_structure(pdb_file_list[j].c_str(),0,1);
                pdb_chain_list[j].second.residues=
                    tmp_model.chains[0].residues;
                tmp_model.chains.clear();
            }
        }
        max_clust_size=fast_clustering(TMsubclust, pdb_name_list,
            tm_fast_mat, tm_full_mat, pdb_chain_list, TMmin, norm);
        if (max_clust_size>1)
        {
            if (TMsubclust.repr_list.size()>1)
            {
                write_TMclust_result(cluster_filename,TMsubclust,
                    pdb_name_list,TMmin,openmode);
            }

            /* clustering at TMmin + TMstep */
            full_clustering(TMsubclust, pdb_name_list, pdb_file_list,
                pdb_chain_list, tm_fast_mat, tm_full_mat, TMmin, TMmax,
                TMstep, cluster_filename, openmode, norm, MinClustSize);
        }

        /* clean up */
        TMsubclust.repr_list.clear();
        TMsubclust.clust_list.clear();
        if (CacheCoor==-1)
        {
            for (l=0;l<TMclust.clust_list[i].size();l++)
            {
                j=TMclust.clust_list[i][l];
                if (pdb_chain_list[j].second.residues.size())
                    pdb_chain_list[j].second.residues.clear();
            }
        }
    }
    return max_clust_size;
}
#endif
