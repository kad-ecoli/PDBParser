/* header for structure alignment by chi-1 rotamer */
#ifndef ROTSUMalign_HPP
#define ROTSUMalign_HPP 1

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <ctype.h>

#include "PDBParser.hpp"
#include "FilePathParser.hpp"
#include "getRotSeq.hpp"
#include "NWalign.hpp"

using namespace std;

void trace_back_gotoh(string seq1, string seq2,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& aln1, string& aln2);
void trace_back_sw(string seq1, string seq2,
    const vector<vector<int> >& JumpH, const vector<vector<int> >& JumpV,
    const vector<vector<int> >& P, string& aln1, string& aln2);
void init_gotoh_mat(vector<vector<int> >&JumpH, vector<vector<int> >&JumpV,
    vector<vector<int> >& P,vector<vector<int> >& S, 
    vector<vector<int> >& H, vector<vector<int> >& V,
    const int len1, const int len2, const int gapopen,const int gapext,
    const int glocal, const int alt_init);
void find_highest_align_score(
    const vector<vector<int> >& S, vector<vector<int> >& P,
    int &aln_score, const int len1,const int len2);

// in rotsum
const int gapopen_rotsum8=-11;
const int gapext_rotsum8=-1;

const int ROTSUM8[55][55]={
//A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z  0  1  2
{ 4, 0,-2, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0, 0, 0, 0,-1,-1, 0,-1, 0, 0, 0, 0,-1,-1,-1, 0,-1,-1, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2,-3,-3,-1,-1,-2},//A
{ 0,12, 8, 4,-3,-5,-6,-1, 0,-1,-1,-4, 0,-5,-1,-6,-4,-1,-3,-5,-4,-9,-6,-4,-5,-6,-2,-6,-6,-1,-4, 0,-6,-6,-5,-5,-5,-5,-1,-5,-3,-5,-5,-3,-2,-6,-2,-2,-1,-2,-7,-9,-3,-6, 2},//B
{-2, 8,11, 7,-6,-5,-6,-2, 0,-3,-3,-3,-1,-6,-4,-4,-4,-4,-4,-3,-6,-6,-6,-6,-1,-4,-4,-2,-5,-3,-3,-1,-5,-5,-7,-6,-6,-6,-2,-5,-4,-2,-4,-2,-3,-2,-3,-2,-4,-4,-5,-8,-6,-5, 0},//C
{ 0, 4, 7,11,-3,-3,-3,-1,-2,-2,-1,-1, 0,-2,-2,-3,-2, 0,-3, 0,-2,-4,-3,-5,-1, 0,-1,-1, 1,-2,-1,-1,-4,-3,-3,-5,-3,-4,-2,-1,-2,-2,-1,-2,-3,-2,-2, 0,-1,-3,-4,-3,-3,-3,-2},//D
{-1,-3,-6,-3,10, 3, 2, 2, 0, 0,-3,-3, 0, 0, 1,-3,-1,-2,-3,-4, 0,-3,-1,-2,-5,-3,-2,-4,-3, 4, 0, 0,-1,-1, 0,-2, 0, 0,-2, 0, 0,-1,-1, 1, 0,-2,-2,-4,-2, 0,-2,-2, 0,-1, 0},//E
{-1,-5,-5,-3, 3, 9, 1, 0, 1, 0,-2,-2,-3,-1, 0, 0,-1,-2,-2,-3, 0,-2, 0,-2,-2,-3,-3,-2,-3, 0, 3, 0,-1,-1, 0, 0, 0,-1, 0, 0, 0, 1,-1,-1, 0,-2,-3,-4,-3,-2,-1,-5, 0, 0,-3},//F
{ 0,-6,-6,-3, 2, 1, 7, 3, 1, 2, 0,-2,-3, 0, 0,-1, 0,-2,-2,-3, 0,-1, 0,-5,-4,-3,-3,-2,-2, 0, 0, 1,-1, 0, 0,-1, 0,-2,-1, 0, 0,-1, 0,-1,-2,-1,-2,-2,-2,-3,-4,-4,-2,-2,-3},//G
{-1,-1,-2,-1, 2, 0, 3, 8, 4, 4, 0,-3,-4,-1, 1,-2, 0,-2,-2,-2, 1, 0, 0,-3,-3,-3, 0,-2,-2, 0,-1, 0,-1,-1, 1, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1, 0,-1,-2,-1, 0,-3,-2, 2,-2,-4},//H
{ 0, 0, 0,-2, 0, 1, 1, 4, 6, 4,-2,-1,-4,-2, 0, 0, 0,-2,-3,-2, 0, 1, 0,-3,-2,-3,-2,-1,-1, 0, 0, 0,-1,-1, 0, 1, 0,-1, 1, 0,-1,-1,-1,-2,-1, 0,-2,-2,-2,-2,-1,-5,-1,-1,-4},//I
{ 0,-1,-3,-2, 0, 0, 2, 4, 4, 5, 0,-2,-2,-1, 0,-1, 0,-2,-2,-3, 0, 0, 0,-2,-3,-2,-1,-2,-1, 0,-1, 0, 0,-1, 0, 0, 2,-1, 0, 0, 0,-1,-1,-1,-2,-1,-2,-3,-1,-3,-1,-3,-4,-2,-3},//J
{-1,-1,-3,-1,-3,-2, 0, 0,-2, 0,12, 1, 0,-1, 0,-4,-1, 1, 0, 0,-3,-5,-5, 0,-2,-2, 1,-1,-1,-1,-2,-2,-4,-4, 0,-4,-3,-3,-4,-3,-1,-2,-3,-3,-1,-4, 0,-1, 0, 5,-3,-1, 8, 1,-1},//K
{-1,-4,-3,-1,-3,-2,-2,-3,-1,-2, 1, 9, 0,-3,-2, 0, 0, 1, 1, 0,-2,-2,-3,-2, 2, 0,-2, 2, 0,-3,-1,-3,-2,-2,-1,-2,-2,-3,-1,-2,-2,-1,-3,-3,-2,-2,-1, 0,-1, 0, 3,-2, 0, 7,-2},//L
{ 0, 0,-1, 0, 0,-3,-3,-4,-4,-2, 0, 0, 7,-2,-2,-4, 0, 0, 0, 0,-3,-5,-3,-2, 0, 1,-1,-1, 1,-3,-4,-2,-1,-1,-3,-5,-3,-5,-4,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0,-1,-3, 3,-2,-2, 3},//M
{ 0,-5,-6,-2, 0,-1, 0,-1,-2,-1,-1,-3,-2, 5, 0,-2, 0,-1,-1,-2, 0,-2,-1,-1,-3,-3,-1,-3,-2, 0, 0, 0,-1,-2, 1,-3, 0, 0,-1,-1, 0, 0, 0,-1,-1,-2,-2,-2,-2,-1,-2,-4, 0,-2,-3},//N
{ 0,-1,-4,-2, 1, 0, 0, 1, 0, 0, 0,-2,-2, 0,14, 3, 2, 0,-1,-3, 0,-2, 0,-1,-4,-2, 0,-1, 0, 3, 0, 0, 0, 0, 3,-1, 0, 2, 0, 0, 0,-1,-2, 0, 0,-2,-1,-3,-2, 3, 0,-3, 5, 0, 0},//O
{ 0,-6,-4,-3,-3, 0,-1,-2, 0,-1,-4, 0,-4,-2, 3,11, 2,-2,-3,-2, 0, 0, 0,-4,-2,-3,-3,-1,-1,-1, 1,-1,-1, 0, 0, 2, 0,-2, 0, 0,-2,-1,-3,-2,-1,-1,-4,-2,-3, 0, 0,-5,-1, 1,-2},//P
{-1,-4,-4,-2,-1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 2, 2, 9, 0,-3,-2,-1,-2, 0,-4,-4,-1,-2,-2, 0, 0, 0, 2, 0,-1, 0, 0, 1,-2,-1, 0,-1, 0, 0, 0, 0, 0,-3,-3,-1,-1,-2,-1,-1,-1, 1},//Q
{-1,-1,-4, 0,-2,-2,-2,-2,-2,-2, 1, 1, 0,-1, 0,-2, 0, 9, 3, 2,-1,-3,-2, 1, 0, 2, 3, 0, 0, 0,-1,-2,-2,-2, 0,-2,-1,-2,-1,-1,-2,-1,-2, 0, 0,-2, 1, 0, 4, 1,-2,-2, 0, 0,-1},//R
{ 0,-3,-4,-3,-3,-2,-2,-2,-3,-2, 0, 1, 0,-1,-1,-3,-3, 3,10, 1,-3,-3,-2, 0, 1, 0, 1, 1, 1,-3,-2,-3, 0,-2,-4,-3,-3,-3,-2,-2,-2, 0,-2,-3, 0,-3, 3, 0, 1, 0,-1, 0, 0, 0,-1},//S
{-1,-5,-3, 0,-4,-3,-3,-2,-2,-3, 0, 0, 0,-2,-3,-2,-2, 2, 1, 6,-3,-2,-2,-1, 0, 2, 2, 1, 2,-3,-3,-2,-4,-4,-4,-2,-2,-4,-1,-2,-3,-2,-3,-2,-1, 0, 0, 2, 0,-2,-2,-2,-2,-1,-1},//T
{ 0,-4,-6,-2, 0, 0, 0, 1, 0, 0,-3,-2,-3, 0, 0, 0,-1,-1,-3,-3, 7, 4, 3,-2,-3,-3,-1,-2,-2, 1, 1, 0, 0, 0, 1, 0, 0, 3, 2, 1, 0,-1, 0, 0, 0,-1,-1,-2,-1,-1,-3,-3,-1,-2,-4},//U
{ 0,-9,-6,-4,-3,-2,-1, 0, 1, 0,-5,-2,-5,-2,-2, 0,-2,-3,-3,-2, 4, 6, 3,-4,-1,-3,-3,-1,-2,-1, 1, 0,-1, 0, 0, 1, 0, 3, 2, 1,-1, 0, 0, 0, 0, 0,-1,-1,-1,-4,-3,-4,-5,-2,-5},//V
{ 0,-6,-6,-3,-1, 0, 0, 0, 0, 0,-5,-3,-3,-1, 0, 0, 0,-2,-2,-2, 3, 3, 5,-4,-3,-1,-1,-2, 0, 0, 0, 1,-1, 0, 0, 0, 2, 0, 0, 2, 0,-1, 0, 0, 0, 0,-2,-2,-2,-4,-3,-1,-5,-3,-2},//W
{ 0,-4,-6,-5,-2,-2,-5,-3,-3,-2, 0,-2,-2,-1,-1,-4,-4, 1, 0,-1,-2,-4,-4, 9, 4, 2, 6, 5,-1,-2,-4,-4,-4,-3,-2,-2,-2,-2,-2,-3,-2,-1, 1,-3,-2,-2, 1, 1, 2, 0,-2,-5, 1,-2,-2},//X
{-1,-5,-1,-1,-5,-2,-4,-3,-2,-3,-2, 2, 0,-3,-4,-2,-4, 0, 1, 0,-3,-1,-3, 4, 6, 1, 2, 4, 0, 0,-2,-3,-3,-3,-3,-1,-2,-3, 0,-2,-3,-2,-1,-3,-2,-1, 3, 1, 0,-2, 2,-3,-3, 0,-2},//Y
{-1,-6,-4, 0,-3,-3,-3,-3,-3,-2,-2, 0, 1,-3,-2,-3,-1, 2, 0, 2,-3,-3,-1, 2, 1, 5, 1, 0, 2,-2,-2,-2,-2, 0,-2,-2,-1,-1,-1, 0,-3,-2,-1, 0,-1,-1, 0, 0, 1,-2, 0,-2,-3, 0,-1},//Z
{-1,-2,-4,-1,-2,-3,-3, 0,-2,-1, 1,-2,-1,-1, 0,-3,-2, 3, 1, 2,-1,-3,-1, 6, 2, 1,11, 6, 4, 0,-2,-3,-2,-2, 0,-2,-1,-1,-2,-1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 2,-3,-3, 3,-1, 0},//a
{ 0,-6,-2,-1,-4,-2,-2,-2,-1,-2,-1, 2,-1,-3,-1,-1,-2, 0, 1, 1,-2,-1,-2, 5, 4, 0, 6, 9, 5,-2,-2,-3,-3,-3,-3,-1,-2,-3, 0,-2,-3,-2,-3,-2,-1, 0, 0, 0, 0, 0, 0,-2,-2, 0,-2},//b
{-1,-6,-5, 1,-3,-3,-2,-2,-1,-1,-1, 0, 1,-2, 0,-1, 0, 0, 1, 2,-2,-2, 0,-1, 0, 2, 4, 5, 8,-2,-2,-2,-2,-3,-1,-2, 0,-3,-2,-1,-3,-3,-2, 0,-1,-1,-1, 0, 0,-2,-1,-1,-3,-1,-1},//c
{-1,-1,-3,-2, 4, 0, 0, 0, 0, 0,-1,-3,-3, 0, 3,-1, 0, 0,-3,-3, 1,-1, 0,-2, 0,-2, 0,-2,-2, 9, 3, 2,-1, 0, 0,-2,-1, 0,-1, 0, 2, 0, 0, 2, 0,-1,-1,-2,-1, 0,-3,-4, 1,-1,-1},//d
{ 0,-4,-3,-1, 0, 3, 0,-1, 0,-1,-2,-1,-4, 0, 0, 1, 0,-1,-2,-3, 1, 1, 0,-4,-2,-2,-2,-2,-2, 3, 8, 3,-2,-1, 0, 0, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0,-3,-3,-2,-2,-2,-3,-1, 0,-1},//e
{-1, 0,-1,-1, 0, 0, 1, 0, 0, 0,-2,-3,-2, 0, 0,-1, 2,-2,-3,-2, 0, 0, 1,-4,-3,-2,-3,-3,-2, 2, 3, 6,-3,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-3,-2,-2,-3,-5, 1,-3,-2, 0},//f
{ 0,-6,-5,-4,-1,-1,-1,-1,-1, 0,-4,-2,-1,-1, 0,-1, 0,-2, 0,-4, 0,-1,-1,-4,-3,-2,-2,-3,-2,-1,-2,-3, 7, 5,-2,-2,-2,-2, 0,-1,-1,-1,-1,-1,-1,-2,-2,-2,-1,-2,-3,-5,-3,-2,-2},//g
{ 0,-6,-5,-3,-1,-1, 0,-1,-1,-1,-4,-2,-1,-2, 0, 0,-1,-2,-2,-4, 0, 0, 0,-3,-3, 0,-2,-3,-3, 0,-1,-2, 5, 7, 0,-1,-1,-1,-1,-1, 0,-1,-1, 0,-1,-2,-2,-2,-1,-1,-3,-5,-3,-3,-3},//h
{ 0,-5,-7,-3, 0, 0, 0, 1, 0, 0, 0,-1,-3, 1, 3, 0, 0, 0,-4,-4, 1, 0, 0,-2,-3,-2, 0,-3,-1, 0, 0, 0,-2, 0, 8, 3, 4, 1, 1,-1, 0, 0,-1, 0, 0, 0,-2,-2,-2, 2,-4,-4, 0,-3,-5},//i
{ 0,-5,-6,-5,-2, 0,-1, 0, 1, 0,-4,-2,-5,-3,-1, 2, 0,-2,-3,-2, 0, 1, 0,-2,-1,-2,-2,-1,-2,-2, 0, 0,-2,-1, 3, 8, 3, 3, 1,-1,-1,-2,-2,-2,-1, 0,-3,-2,-3,-3,-3,-5,-4,-2,-4},//j
{-1,-5,-6,-3, 0, 0, 0, 0, 0, 2,-3,-2,-3, 0, 0, 0, 1,-1,-3,-2, 0, 0, 2,-2,-2,-1,-1,-2, 0,-1, 0, 0,-2,-1, 4, 3, 6, 0, 0, 0, 0,-1, 0, 0,-1, 0,-2,-2,-2,-3,-3,-3,-3,-1,-3},//k
{ 0,-5,-6,-4, 0,-1,-2, 0,-1,-1,-3,-3,-5, 0, 2,-2,-2,-2,-3,-4, 3, 3, 0,-2,-3,-1,-1,-3,-3, 0,-1,-1,-2,-1, 1, 3, 0, 9, 4, 4, 0,-1,-2,-1,-1,-3,-2,-4,-2, 0,-2,-4, 1,-2,-4},//l
{ 0,-1,-2,-2,-2, 0,-1, 0, 1, 0,-4,-1,-4,-1, 0, 0,-1,-1,-2,-1, 2, 2, 0,-2, 0,-1,-2, 0,-2,-1, 0, 0, 0,-1, 1, 1, 0, 4, 7, 3,-1, 0, 0,-2, 0, 0,-1,-2,-1, 0, 0,-4,-1, 0,-3},//m
{ 0,-5,-5,-1, 0, 0, 0, 0, 0, 0,-3,-2,-2,-1, 0, 0, 0,-1,-2,-2, 1, 1, 2,-3,-2, 0,-1,-2,-1, 0, 0, 0,-1,-1,-1,-1, 0, 4, 3, 6,-1, 0, 0, 0,-1,-1,-3,-3,-2,-1,-1,-2,-3,-1,-1},//n
{ 0,-3,-4,-2, 0, 0, 0, 0,-1, 0,-1,-2,-2, 0, 0,-2,-1,-2,-2,-3, 0,-1, 0,-2,-3,-3,-1,-3,-3, 2, 0, 0,-1, 0, 0,-1, 0, 0,-1,-1, 4, 2, 2, 2, 0, 0,-1,-2, 0, 0,-3,-4,-1,-1,-2},//o
{ 0,-5,-2,-2,-1, 1,-1,-1,-1,-1,-2,-1,-2, 0,-1,-1, 0,-1, 0,-2,-1, 0,-1,-1,-2,-2, 0,-2,-3, 0, 1, 0,-1,-1, 0,-2,-1,-1, 0, 0, 2, 6, 2, 0, 1, 0,-2,-2,-1,-3,-3,-3,-1,-1,-1},//p
{ 0,-5,-4,-1,-1,-1, 0,-1,-1,-1,-3,-3,-1, 0,-2,-3, 0,-2,-2,-3, 0, 0, 0, 1,-1,-1, 0,-3,-2, 0, 0, 0,-1,-1,-1,-2, 0,-2, 0, 0, 2, 2, 4, 0, 0, 1,-2,-2,-2,-2,-4,-2,-1,-2,-1},//q
{-1,-3,-2,-2, 1,-1,-1,-1,-2,-1,-3,-3,-1,-1, 0,-2, 0, 0,-3,-2, 0, 0, 0,-3,-3, 0,-1,-2, 0, 2, 0, 0,-1, 0, 0,-2, 0,-1,-2, 0, 2, 0, 0, 5, 3, 2, 0,-1, 0, 0,-3,-3,-2,-2,-1},//r
{ 0,-2,-3,-3, 0, 0,-2,-1,-1,-2,-1,-2,-1,-1, 0,-1, 0, 0, 0,-1, 0, 0, 0,-2,-2,-1, 0,-1,-1, 0, 0, 0,-1,-1, 0,-1,-1,-1, 0,-1, 0, 1, 0, 3, 7, 3, 0,-1, 0,-1,-3,-3,-1,-2,-1},//s
{ 0,-6,-2,-2,-2,-2,-1, 0, 0,-1,-4,-2,-1,-2,-2,-1, 0,-2,-3, 0,-1, 0, 0,-2,-1,-1, 0, 0,-1,-1, 0, 0,-2,-2, 0, 0, 0,-3, 0,-1, 0, 0, 1, 2, 3, 5, 0, 0,-1,-3,-3,-2,-2,-2,-1},//t
{ 0,-2,-3,-2,-2,-3,-2,-1,-2,-2, 0,-1, 0,-2,-1,-4,-3, 1, 3, 0,-1,-1,-2, 1, 3, 0, 0, 0,-1,-1,-3,-3,-2,-2,-2,-3,-2,-2,-1,-3,-1,-2,-2, 0, 0, 0, 7, 2, 4, 0,-2,-2, 0,-2,-1},//u
{ 0,-2,-2, 0,-4,-4,-2,-2,-2,-3,-1, 0, 0,-2,-3,-2,-3, 0, 0, 2,-2,-1,-2, 1, 1, 0, 0, 0, 0,-2,-3,-2,-2,-2,-2,-2,-2,-4,-2,-3,-2,-2,-2,-1,-1, 0, 2, 4, 2,-3,-3,-1,-3,-2, 0},//v
{ 0,-1,-4,-1,-2,-3,-2,-1,-2,-1, 0,-1, 0,-2,-2,-3,-1, 4, 1, 0,-1,-1,-2, 2, 0, 1, 0, 0, 0,-1,-2,-2,-1,-1,-2,-3,-2,-2,-1,-2, 0,-1,-2, 0, 0,-1, 4, 2, 7,-1,-3,-2,-2,-2,-1},//w
{-2,-2,-4,-3, 0,-2,-3, 0,-2,-3, 5, 0,-1,-1, 3, 0,-1, 1, 0,-2,-1,-4,-4, 0,-2,-2, 2, 0,-2, 0,-2,-3,-2,-1, 2,-3,-3, 0, 0,-1, 0,-3,-2, 0,-1,-3, 0,-3,-1,16, 4, 0, 7, 1,-2},//x
{-3,-7,-5,-4,-2,-1,-4,-3,-1,-1,-3, 3,-3,-2, 0, 0,-2,-2,-1,-2,-3,-3,-3,-2, 2, 0,-3, 0,-1,-3,-2,-5,-3,-3,-4,-3,-3,-2, 0,-1,-3,-3,-4,-3,-3,-3,-2,-3,-3, 4,12, 0, 1, 6,-4},//y
{-3,-9,-8,-3,-2,-5,-4,-2,-5,-3,-1,-2, 3,-4,-3,-5,-1,-2, 0,-2,-3,-4,-1,-5,-3,-2,-3,-2,-1,-4,-3, 1,-5,-5,-4,-5,-3,-4,-4,-2,-4,-3,-2,-3,-3,-2,-2,-1,-2, 0, 0,11,-2,-2, 0},//z
{-1,-3,-6,-3, 0, 0,-2, 2,-1,-4, 8, 0,-2, 0, 5,-1,-1, 0, 0,-2,-1,-5,-5, 1,-3,-3, 3,-2,-3, 1,-1,-3,-3,-3, 0,-4,-3, 1,-1,-3,-1,-1,-1,-2,-1,-2, 0,-3,-2, 7, 1,-2,11, 4, 0},//0
{-1,-6,-5,-3,-1, 0,-2,-2,-1,-2, 1, 7,-2,-2, 0, 1,-1, 0, 0,-1,-2,-2,-3,-2, 0, 0,-1, 0,-1,-1, 0,-2,-2,-3,-3,-2,-1,-2, 0,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2, 1, 6,-2, 4, 9, 0},//1
{-2, 2, 0,-2, 0,-3,-3,-4,-4,-3,-1,-2, 3,-3, 0,-2, 1,-1,-1,-1,-4,-5,-2,-2,-2,-1, 0,-2,-1,-1,-1, 0,-2,-3,-5,-4,-3,-4,-3,-1,-2,-1,-1,-1,-1,-1,-1, 0,-1,-2,-4, 0, 0, 0, 8},//2
};


const string RotSeq_list=
//ACCCDDDEEEFFFGHHHIIIKKKLLLMMMNNNPPQQQRRRSSSTTTVVVWWWYYY amino acid code
//1123123123123112312312312312312312123123123123123123123 chi-1 rotamer index
 "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz012"; // jarret rotamer code

/* convert ROTSUM code to int */
inline int RotSeq2int(char aa)
{
    for (int i=0;i<RotSeq_list.length();i++)
        if (RotSeq_list[i]==aa) return i;
    //if (aa!=toupper(aa)) return RotSeq2int(toupper(aa));
    return RotSeq_list.length();
}

vector<int> RotSeq2int(const string sequence)
{
    vector<int> seq2int;
    for (int r=0;r<sequence.length();r++)
        seq2int.push_back(RotSeq2int(sequence[r]));
    return seq2int;
}

/* calculate dynamic programming matrix using gotoh algorithm. 
 * overwriting NWalign in NWalign.hpp when scoring matrix is int 55x55.
 *
 * S     - cumulative scorefor each cell
 * P     - string representation for path
 *         0 :   uninitialized, for gaps at N- & C- termini when glocal>0
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seq1)+1]*[len(seq2)+1]
 *
 * global - global or local alignment
 *         0 : global alignment (Needleman-Wunsch dynamic programming)
 *         1 : glocal-query alignment
 *         2 : glocal-both alignment
 *         3 : local alignment (Smith-Waterman dynamic programming)
 *
 * alt_init - whether to adopt alternative matrix initialization
 *         1 : use wei zheng's matrix initialization
 *         0 : use yang zhang's matrix initialization, does NOT work
 *             for glocal alignment
 */
int calculate_score_gotoh(
    const vector<int>& seq2int1, const vector<int>& seq2int2,
    vector<vector<int> >& JumpH, vector<vector<int> >& JumpV,
    vector<vector<int> >& P,const int ScoringMatrix[55][55],
    const int gapopen,const int gapext,const int glocal=0,
    const int alt_init=1)
{
    int len1=seq2int1.size();
    int len2=seq2int2.size();

    vector<int> temp_int(len2+1,0);
    vector<vector<int> > S(len1+1,temp_int);
    // penalty score for horizontal long gap
    vector<vector<int> > H(len1+1,temp_int);
    // penalty score for vertical long gap
    vector<vector<int> > V(len1+1,temp_int);
    
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i,j;
    init_gotoh_mat(JumpH, JumpV, P, S, H, V, len1, len2,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score,left_score,up_score;
    for (i=1;i<len1+1;i++)
    {
        for (j=1;j<len2+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<len1 || glocal>=3)
            {
                H[i][j]=MAX(S[i][j-1]+gapopen,H[i][j-1]+gapext);
                JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            }
            else
            {
                H[i][j]=MAX(S[i][j-1],H[i][j-1]);
                JumpH[i][j]=(H[i][j]==H[i][j-1])?(JumpH[i][j-1]+1):1;
            }
            // penalty of consective insertion
            if (glocal<2 || j<len2 || glocal>=3)
            {
                V[i][j]=MAX(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=MAX(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]; // match-mismatch '\'
            if (seq2int1[i-1]<55 && seq2int2[j-1]<55)
                diag_score+=ScoringMatrix[seq2int1[i-1]][seq2int2[j-1]];
            left_score=H[i][j];     // deletion       '-'
            up_score  =V[i][j];     // insertion      '|'

            if (diag_score>=left_score && diag_score>=up_score)
            {
                S[i][j]=diag_score;
                P[i][j]+=1;
            }
            if (up_score>=diag_score && up_score>=left_score)
            {
                S[i][j]=up_score;
                P[i][j]+=2;
            }
            if (left_score>=diag_score && left_score>=up_score)
            {
                S[i][j]=left_score;
                P[i][j]+=4;
            }
            if (glocal>=3 && S[i][j]<0)
            {
                S[i][j]=0;
                P[i][j]=0;
                H[i][j]=0;
                V[i][j]=0;
                JumpH[i][j]=0;
                JumpV[i][j]=0;
            }
        }
    }
    int aln_score=S[len1][len2];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<len1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<len2;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(S,P,aln_score,len1,len2);

    // release memory
    S.clear();
    H.clear();
    V.clear();
    return aln_score; // final alignment score
}

#endif
