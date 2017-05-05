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

const int gapopen_3dblast=-8;
const int gapext_3dblast=-2;

const int BLOSUM62_3dblast[24][24]={
// A   B   C   D   E   F  G   H  I  K  L  M  N  P  Q  R  S  T  V  W  X   Y  Z U
{  5,  2,  2,  2,-12,-12,-1, -9,-2,-8, 0,-3,-7,-7,-3,-5,-5,-7,-6,-4,-6,  3,-4,0},//A
{  2,  5,  2,  2,-12,-10, 1,-10,-2,-7,-2,-2,-7,-6,-3,-5,-5,-6,-5,-4,-6,  2,-4,0},//B
{  2,  2,  5,  1,-11, -9,-1, -9, 1,-8,-1,-3,-7,-6,-3,-5,-5,-7,-6,-5,-6,  3,-4,0},//C
{  2,  2,  1,  5,-10, -9, 1, -9, 0,-6, 1,-1,-5,-5,-2,-4,-4,-5,-4,-1,-4,  2,-3,0},//D
{-12,-12,-11,-10,  6,  1,-8,  2,-9,-2,-8,-6,-1,-4,-7,-6,-8,-4,-4,-6,-3,-15,-3,0},//E
{-12,-10, -9, -9,  1,  6,-6,  0,-7, 1,-7,-4,-1,-3,-5,-4,-6,-3,-4,-5,-2,-10,-2,0},//F
{ -1,  1, -1,  1, -8, -6, 7, -5, 0,-4,-1, 2,-4,-3, 1,-2,-3,-3,-1,-1,-2, -1,-2,0},//G
{ -9,-10, -9, -9,  2,  0,-5,  6,-6,-1,-6,-4, 2,-2,-4,-2,-6,-3,-3,-4, 0,-10,-2,0},//H
{ -2, -2,  1,  0, -9, -7, 0, -6, 9,-5, 3,-1,-3,-4,-1,-2,-2,-4,-3, 2,-3, -2,-2,0},//I
{ -8, -7, -8, -6, -2,  1,-4, -1,-5, 6,-6,-3, 1,-3,-4,-4,-4,-1,-2,-4,-1, -8, 0,0},//K
{  0, -2, -1,  1, -8, -7,-1, -6, 3,-6, 7,-2,-5,-4,-1,-1,-1,-3,-2, 3,-4, -1,-1,0},//L
{ -3, -2, -3, -1, -6, -4, 2, -4,-1,-3,-2, 7,-3,-3,-1,-2,-4,-2, 2,-1, 2, -3,-2,0},//M
{ -7, -7, -7, -5, -1, -1,-4,  2,-3, 1,-5,-3, 6, 1,-2, 0,-3, 1,-1,-3, 0, -8, 0,0},//N
{ -7, -6, -6, -5, -4, -3,-3, -2,-4,-3,-4,-3, 1, 7,-2, 1, 0, 1,-2,-2,-2, -7,-1,0},//P
{ -3, -3, -3, -2, -7, -5, 1, -4,-1,-4,-1,-1,-2,-2, 6, 3,-2,-2,-3,-1,-1, -3,-2,0},//Q
{ -5, -5, -5, -4, -6, -4,-2, -2,-2,-4,-1,-2, 0, 1, 3, 8,-2,-1,-2,-2, 1, -5,-2,0},//R
{ -5, -5, -5, -4, -8, -6,-3, -6,-2,-4,-1,-4,-3, 0,-2,-2, 8, 0,-1, 2,-3, -5,-2,0},//S
{ -7, -6, -7, -5, -4, -3,-3, -3,-4,-1,-3,-2, 1, 1,-2,-1, 0, 6, 0,-1,-1, -7,-2,0},//T
{ -6, -5, -6, -4, -4, -4,-1, -3,-3,-2,-2, 2,-1,-2,-3,-2,-1, 0, 8, 2, 1, -7,-1,0},//V
{ -4, -4, -5, -1, -6, -5,-1, -4, 2,-4, 3,-1,-3,-2,-1,-2, 2,-1, 2,11,-2, -6,-2,0},//W
{ -6, -6, -6, -4, -3, -2,-2,  0,-3,-1,-4, 2, 0,-2,-1, 1,-3,-1, 1,-2, 7, -7, 0,0},//X
{  3,  2,  3,  2,-15,-10,-1,-10,-2,-8,-1,-3,-8,-7,-3,-5,-5,-7,-7,-6,-7,  5,-4,0},//Y
{ -4, -4, -4, -3, -3, -2,-2, -2,-2, 0,-1,-2, 0,-1,-2,-2,-2,-2,-1,-2, 0, -4, 9,0},//Z
{  0,  0,  0,  0,  0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,0},//U
};

const string ThreeDblast_list="ABCDETKVNFGHILMQSYRPWZXU"; // U is chain terminal

/* convert ThreeDblast structure alphabet to int */
inline int ThreeDblast2int(char aa)
{
    for (int i=0;i<ThreeDblast_list.length();i++)
        if (ThreeDblast_list[i]==aa) return i;
    if (aa!=toupper(aa)) return aa2int(toupper(aa));
    return ThreeDblast_list.length();
}

vector<int> ThreeDblast2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(ThreeDblast2int(sequence[r]));
    return seq2int;
}

/* read multiple-chain PDB and extract the 3d-blast structure alphabet */
int read_pdb_as_3dblast(const char *filename,vector<string>& name_list,
    vector<string>& seq_list, vector<vector<int> >& seq2int_list)
{
    int atomic_detail=0; // only read CA
    int allowX=1;        // only allow ATOM and MSE

    string PDBid=basename_no_ext(filename);
    ModelUnit pdb_entry=read_pdb_structure(filename,atomic_detail,allowX);

    int seq_num=pdb_entry.chains.size();
    string ThreeDblast;
    vector<int> seq2int;
    int c,r;
    for (c=0;c<seq_num;c++)
    {
        ThreeDblast=pdb2ThreeDblast(pdb_entry.chains[c]);
        seq_list.push_back(ThreeDblast);
        seq2int_list.push_back(aa2int(ThreeDblast));

        ThreeDblast.clear();
        seq2int.clear();
        name_list.push_back(PDBid+':'+pdb_entry.chains[c].chainID_full);
    }
    pdb_entry.chains.clear();
    return seq_num;
}
