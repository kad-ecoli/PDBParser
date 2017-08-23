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
    vector<int> unclust_list; // index of unclustered items
    vector<int> repr_list;    // index of cluster representative
    vector<vector<int> > clust_list; // each row is one cluster
                              // each column is one member in a cluster
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

    int p,s;
    bool pdb_file_found=false;
    string pdb_file_fullname;
    while(fp.good())
    {
        getline(fp,line);
        if (line.length()==0 || line[0]=='#') continue;

        bool pdb_file_found=false;
        for (p=0;p<prefix_list.size()&&pdb_file_found==false;p++)
        {
            for (s=0;s<suffix_list.size()&&pdb_file_found==false;s++)
            {
                pdb_file_fullname=pdb_folder_fullname+
                    prefix_list[p]+line+suffix_list[s];
                if (isfile(pdb_file_fullname.c_str()))
                {
                    pdb_file_list.push_back(pdb_file_fullname);
                    pdb_file_found=true;
                    break;
                }
            }
        }

        if (pdb_file_found) pdb_name_list.push_back(line);
        else cerr<<"ERROR! Cannot locate PDB file for "<<line<<endl;
    }
    fp.close();
    return pdb_name_list.size();
}

/* function for sorting vector<pair<int,ChainUnit> > */ 
bool cf_pdb_chain_list(pair<int,ChainUnit> chain1,pair<int,ChainUnit> chain2)
{
    return chain1.first>chain2.first;
}

bool cf_int_str_pair_list(pair<int,string> chain1,pair<int,string> chain2)
{
    return chain1.first>chain2.first;
}

void batch_sarst_rmsd(const vector<string> &pdb_name_list,
    vector<pair<int,ChainUnit> >&pdb_chain_list,
    vector<vector<int> >&sarst2int_list,
    vector<vector<double> >&tm_fast_mat, vector<vector<double> >&tm_full_mat,
    double tmscore_cutoff=0.5, int norm=0)
{
    vector<double> tmp_array(3,0);
    string aln_i,aln_j;
    int aln_len;
    int pdb_entry_num=pdb_chain_list.size();

    vector<vector<double> > xyz_list_i,xyz_list_j;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t

    double rmsd,tmscore_i,tmscore_j;

    int i,j;
    //for (i=0;i<pdb_entry_num-1;i++)
    for (i=0;i<pdb_entry_num;i++)
    {
        cout<<"    aligning "<<pdb_name_list[i]<<endl;
        //for (j=i+1;j<pdb_entry_num;j++)
        for (j=0;j<i;j++)
        {
            NWalign(pdb_chain_list[i].second.sarst, pdb_chain_list[j].
                second.sarst, sarst2int_list[i], sarst2int_list[j],
                aln_i, aln_j, BLOSUM62_sarst, gapopen_sarst,gapext_sarst,0);

            aln2coor(aln_i, aln_j, pdb_chain_list[i].second,
                pdb_chain_list[j].second, xyz_list_i, xyz_list_j,0);
            aln_len=xyz_list_i.size();

            if (aln_len==0) continue;

            /* kabsch */
            RotateCoor(xyz_list_i,xyz_list_j, RotMatix, TranVect);

            /* change coordinate */
            vector<vector<double> > super_xyz_list_i(aln_len,tmp_array);
            for (int r=0; r<aln_len; r++)
	            ChangeCoor(xyz_list_i[r], RotMatix, TranVect, 
                    super_xyz_list_i[r]);

            /* TM-score */
            tmscore_i=calTMscore(super_xyz_list_i,xyz_list_j,
                pdb_chain_list[i].first);
            tmscore_j=calTMscore(super_xyz_list_i,xyz_list_j,
                pdb_chain_list[j].first);
                
            /* clean up */
            aln_i.clear();
            aln_j.clear();
            RotMatix.clear();
            TranVect.clear();
            super_xyz_list_i.clear();
            xyz_list_i.clear();
            xyz_list_j.clear();

            /* perform TMalign on well aligned pairs */
            if (tmscore_j>=tmscore_cutoff*0.6||(norm==1 && 
                tmscore_i>=tmscore_cutoff*0.6))
            {
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, rmsd, tmscore_i, tmscore_j,0);
                tm_fast_mat[i][j]=tm_full_mat[i][j]=tmscore_i;
                tm_fast_mat[j][i]=tm_full_mat[j][i]=tmscore_j;
                cout<<pdb_name_list[i]<<'\t'<<pdb_name_list[j]<<'\t'
                    <<setiosflags(ios::fixed)<<setprecision(4)
                    <<tmscore_i<<'\t'<<tmscore_j<<endl;
            }
        }
    }
}

void initialize_TMclust(TMclustUnit &TMclust,const int pdb_entry_num)
{
    TMclust.repr_list.clear();
    TMclust.clust_list.clear();
    TMclust.unclust_list.clear();

    vector<int> tmp_array(1,0);
    TMclust.repr_list.push_back(0); // first chain is first cluster
    TMclust.clust_list.push_back(tmp_array);
    for (int i=pdb_entry_num-1;i>0;i--) // vector cannot pop_front. so reverse
        TMclust.unclust_list.push_back(i); // the order of member list
}

/* initialize TMclustUnit with members list */
void initialize_TMclust(TMclustUnit &TMclust, const vector<int> &member_list)
{
    TMclust.repr_list.clear();
    TMclust.clust_list.clear();
    TMclust.unclust_list.clear();
    if (member_list.size()==0) return;

    vector<int> tmp_array(1,member_list[0]);
    TMclust.repr_list.push_back(member_list[0]); // first chain is first cluster
    TMclust.clust_list.push_back(tmp_array);
    for (int i=member_list.size()-1;i>0;i--)
        TMclust.unclust_list.push_back(member_list[i]);
}

/* add new_member to cluster represented by clust_repr */
int add_chain_to_clust(TMclustUnit &TMclust, const int new_member,
    const int clust_repr)
{
    int j;
    for (j=0;j<TMclust.repr_list.size();j++)
    {
        if (TMclust.repr_list[j]==clust_repr)
        {
            TMclust.clust_list[j].push_back(new_member);
            break;
        }
    }
    if (j==TMclust.repr_list.size())
    {
        cout<<"FATAL ERROR! No such representative "<<clust_repr<<endl;
        exit(0);
    }
    return j; // index of clust_repr in repr_list
}

int fast_clustering(TMclustUnit &TMclust,const vector<string>&pdb_name_list,
    vector<vector<double> >&tm_fast_mat, vector<vector<double> >&tm_full_mat,
    vector<pair<int,ChainUnit> >&pdb_chain_list, double tmscore_cutoff=0.5,
    int f=8, int norm=0)
{
    int i,j,k;
    int add_to_clust=-1;
    int max_clust_size=0;

    string aln_i,aln_j;
    double rmsd,tmscore_i,tmscore_j;
    while(TMclust.unclust_list.size())
    {
        i=TMclust.unclust_list.back(); // try to add protein i to clust_list
        cout<<"assigning "<<pdb_name_list[i];
        TMclust.unclust_list.pop_back();
        add_to_clust=-1; // cluster representative whose cluster i belongs
        for (k=0;k<TMclust.repr_list.size();k++)
        {
            j=TMclust.repr_list[k];
            if (MIN(tm_full_mat[i][j],tm_full_mat[j][i])>=tmscore_cutoff
                || (norm==1 && 
                MAX(tm_full_mat[i][j],tm_full_mat[j][i])>=tmscore_cutoff))
            {
                add_to_clust=j;
                break;
            }
        }
        if (add_to_clust>=0) 
        {
            cout<<" to "<<pdb_name_list[j]<<endl;
            add_chain_to_clust(TMclust,i,j);
            continue;
        }
        for (k=0;k<TMclust.repr_list.size();k++)
        {
            j=TMclust.repr_list[k];
            if (tm_fast_mat[i][j]!=0) continue;
        
            rmsd=tmscore_i=tmscore_j=0;
            TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                pdb_chain_list[j].second, rmsd, tmscore_i, tmscore_j, f);
            tm_fast_mat[i][j]=tmscore_i;
            tm_fast_mat[j][i]=tmscore_j;
            if (f>0     &&(MIN(tmscore_i,tmscore_j)>=tmscore_cutoff*.9 ||
               (norm==1 && MAX(tmscore_i,tmscore_j)>=tmscore_cutoff*.9)))
            {
                TMalign(aln_i, aln_j,pdb_chain_list[i].second,
                    pdb_chain_list[j].second, rmsd, tmscore_i, tmscore_j, 0);
                tm_full_mat[i][j]=tmscore_i;
                tm_full_mat[j][i]=tmscore_j;
            }
            else if (f==0)
            {
                tm_full_mat[i][j]=tmscore_i;
                tm_full_mat[j][i]=tmscore_j;
            }
            if (MIN(tmscore_i,tmscore_j)>=tmscore_cutoff|| (norm==1 && 
                MAX(tmscore_i,tmscore_j)>=tmscore_cutoff))
            {
                add_to_clust=j;
                break;
            }
        }
        if (add_to_clust>=0) 
        {
            cout<<" to "<<pdb_name_list[j]<<endl;
            add_chain_to_clust(TMclust,i,j);
            continue;
        }
        TMclust.repr_list.push_back(i);
        vector<int> tmp_array(1,i);
        TMclust.clust_list.push_back(tmp_array);
        cout<<" as new cluster"<<endl;
    }
    for (i=0;i<TMclust.clust_list.size();i++) 
        max_clust_size=MAX(max_clust_size,TMclust.clust_list[i].size());
    return max_clust_size;
}

void write_TMclust_result(const string filename,
    const TMclustUnit &TMclust, const vector<string> pdb_name_list,
    double tmscore_cutoff=0.5,char openmode='w')
{
    stringstream buf;
    int i,j;

    int N=0; // total number of entries
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

void write_matrix(const string filename, const vector<vector<double> >&tm_mat)
{
    stringstream buf;
    int i,j;
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
    const vector<int>&repr_list, const vector<string>& pdb_name_list,
    vector<pair<int,ChainUnit> >& pdb_chain_list,char openmode='w')
{
    stringstream buf;
    int repr,i,r,res_num;
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

string full_clustering(TMclustUnit &TMclust,const vector<string>&pdb_name_list,
    vector<pair<int,ChainUnit> >&pdb_chain_list,
    vector<vector<double> >&tm_fast_mat, vector<vector<double> >&tm_full_mat,
    double TMmin=0.5, double TMmax=0.8, double TMstep=0.1,
    string cluster_filename="cluster.txt", string ca_xyz_filename="ca.xyz",
    char openmode='a', int norm=0, int MinClustSize=2)
{
    string clust_txt="";
    TMmin+=TMstep;
    if (TMmin>TMmax) return clust_txt;
    if (TMclust.repr_list.size()<=0) return clust_txt;

    cout<<"heuristic clustering for TM-score "<<TMmin<<endl;
    int max_clust_size=0;
    for (int i=0;i<TMclust.repr_list.size();i++)
    {
        if (TMclust.clust_list[i].size()<=MinClustSize) continue;

        /* clustering at TMmin */
        TMclustUnit TMsubclust;
        initialize_TMclust(TMsubclust,TMclust.clust_list[i]);
        max_clust_size=fast_clustering(TMsubclust, pdb_name_list,
            tm_fast_mat, tm_full_mat, pdb_chain_list, TMmin, 0, norm);
        if (max_clust_size<=1) continue;
        if (TMsubclust.repr_list.size()>1)
        {
            write_TMclust_result(cluster_filename,TMsubclust,
                pdb_name_list,TMmin,openmode);
        }

        //vector<int> repr_list;
        //for (int repr=1;repr<TMsubclust.repr_list.size();repr++)
            //repr_list.push_back(TMsubclust.repr_list[repr]);
        //write_TMclust_ca_xyz(ca_xyz_filename, repr_list,
        //    pdb_name_list, pdb_chain_list, openmode);

        /* clustering at TMmin + TMstep */
        full_clustering(TMsubclust, pdb_name_list, pdb_chain_list,
            tm_fast_mat, tm_full_mat, TMmin, TMmax, TMstep,
            cluster_filename, ca_xyz_filename, openmode, norm);

        /* clean up */
        //repr_list.clear();
        TMsubclust.repr_list.clear();
        TMsubclust.clust_list.clear();
    }
    return clust_txt;
}
#endif
