#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <ctype.h>

#include "FilePathParser.hpp"
#include "StructuralAlphabet.hpp"

using namespace std;

/*
// in CPSARST
const int gapopen_sarst=-10;
const int gapext_sarst=-3;
*/

// in SARST
const int gapopen_sarst=-9;
const int gapext_sarst=-2;

const int BLOSUM62_sarst[24][24]={
// A   B   C  D  E  T  K  V  N  F   G   H   I  L  M  Q   S   Y  R  P  W  Z  X  *
{  3,  2,  2, 1, 1, 0,-2,-3,-3,-8,-11,-11,-13,-8,-8,-9,-14, -9,-7,-8,-7,-4, 0,-4},//A
{  2,  2,  2, 1, 1, 1, 0,-1,-2,-6,-12,-10,-10,-7,-7,-6,-10, -8,-5,-6,-4,-6, 0,-4},//B
{  2,  2,  2, 1, 1, 3,-1,-2,-3,-6,-13,-11, -9,-7,-8,-7, -9,-10,-2,-7,-5,-3, 0,-4},//C
{  1,  1,  1, 3, 1, 2, 2,-1,-1,-4, -9, -7, -8,-4,-6,-5, -7, -4, 1,-3,-4,-2, 0,-4},//D
{  1,  1,  1, 1, 3, 1, 2, 3, 1,-5, -7, -6, -7,-4,-4,-4, -7, -2,-1,-5,-3,-1, 0,-4},//E
{  0,  1,  3, 2, 1, 5,-1, 2,-1,-2, -6, -6, -4,-4,-5,-2, -4, -4, 2,-1,-1, 3, 0,-4},//T
{ -2,  0, -1, 2, 2,-1, 4, 1, 3,-3, -6, -6, -5,-3,-3,-2, -5, -2,-2, 0, 0,-1, 0,-4},//K
{ -3, -1, -2,-1, 3, 2, 1, 9, 3,-3, -4, -4, -2,-2,-2, 0,  0,  3, 2,-1, 3, 4, 0,-4},//V
{ -3, -2, -3,-1, 1,-1, 3, 3, 5,-2, -4, -4, -3,-2, 0,-2, -3, -2,-1, 1, 1, 1, 0,-4},//N
{ -8, -6, -6,-4,-5,-2,-3,-3,-2, 5, -1,  1,  0, 3, 0, 3,  0,  2, 0,-2,-2, 1, 0,-4},//F
{-11,-12,-13,-9,-7,-6,-6,-4,-4,-1,  4,  3,  3, 0, 2, 0,  1, -3,-5,-5,-6,-2, 0,-4},//G
{-11,-10,-11,-7,-6,-6,-6,-4,-4, 1,  3,  4,  1, 2, 2, 0, -1, -2,-4,-3,-5,-1, 0,-4},//H
{-13,-10, -9,-8,-7,-4,-5,-2,-3, 0,  3,  1,  4, 0, 1, 2,  4,  0,-1,-4,-7,-2, 0,-4},//I
{ -8, -7, -7,-4,-4,-4,-3,-2,-2, 3,  0,  2,  0, 4, 1, 1, -1,  0, 0,-1,-2, 1, 0,-4},//L
{ -8, -7, -8,-6,-4,-5,-3,-2, 0, 0,  2,  2,  1, 1, 4, 0,  1, -1,-4,-2,-2, 1, 0,-4},//M
{ -9, -6, -7,-5,-4,-2,-2, 0,-2, 3,  0,  0,  2, 1, 0, 6,  1,  3, 1,-3,-3, 1, 0,-4},//Q
{-14,-10, -9,-7,-7,-4,-5, 0,-3, 0,  1, -1,  4,-1, 1, 1,  7,  5, 2,-3,-3, 3, 0,-4},//S
{ -9, -8,-10,-4,-2,-4,-2, 3,-2, 2, -3, -2,  0, 0,-1, 3,  5, 10, 7, 2, 2, 7, 0,-4},//Y
{ -7, -5, -2, 1,-1, 2,-2, 2,-1, 0, -5, -4, -1, 0,-4, 1,  2,  7,11, 3, 0, 7, 0,-4},//R
{ -8, -6, -7,-3,-5,-1, 0,-1, 1,-2, -5, -3, -4,-1,-2,-3, -3,  2, 3, 8, 7, 4, 0,-4},//P
{ -7, -4, -5,-4,-3,-1, 0, 3, 1,-2, -6, -5, -7,-2,-2,-3, -3,  2, 0, 7, 9, 5, 0,-4},//W
{ -4, -6, -3,-2,-1, 3,-1, 4, 1, 1, -2, -1, -2, 1, 1, 1,  3,  7, 7, 4, 5, 6, 0,-4},//Z
{  0,  0,  0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0,-4},//X
{ -4, -4, -4,-4,-4,-4,-4,-4,-4,-4, -4, -4, -4,-4,-4,-4, -4, -4,-4,-4,-4,-4,-4,-4},//*
};

const string sarst_list="ABCDETKVNFGHILMQSYRPWZX*";

/* convert sarst code to int */
inline int sarst2int(char aa)
{
    for (int i=0;i<sarst_list.length();i++) if (sarst_list[i]==aa) return i;
    if (aa!=toupper(aa)) return aa2int(toupper(aa));
    return sarst_list.length();
}

vector<int> sarst2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(sarst2int(sequence[r]));
    return seq2int;
}

/* read multiple-chain PDB and extract the sarst code */
int read_pdb_as_sarst(const char *filename,vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list)
{
    int atomic_detail=1; // only read backbone
    int allowX=1;        // only allow ATOM and MSE

    string PDBid=basename_no_ext(filename);
    ModelUnit pdb_entry=read_pdb_structure(filename,atomic_detail,allowX);

    int seq_num=pdb_entry.chains.size();
    string sarst;
    vector<int> seq2int;
    int c,r;
    for (c=0;c<seq_num;c++)
    {
        sarst=pdb2sarst(pdb_entry.chains[c]);
        seq_list.push_back(sarst);
        seq2int_list.push_back(aa2int(sarst));

        sarst.clear();
        seq2int.clear();
        name_list.push_back(PDBid+':'+pdb_entry.chains[c].chainID_full);
    }
    pdb_entry.chains.clear();
    return seq_num;
}
