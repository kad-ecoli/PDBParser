/* ccealign -- structural alignment by combinatorial extension (CE)
 * 
 * Copyright (c) 2017, Chengxin Zhang.
 * Copyright (c) 2007, Jason Vertrees.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *      
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ccealign_HPP
#define ccealign_HPP 1

#include <math.h>

// include the template numerical toolkit
#include "tnt/tnt.h"

// for reflections
#include "jama_lu.h"

// for the SVD
#include "jama_svd.h"

#include "PDBParser.hpp"

using namespace std;

// CE constants
const float windowSize = 8.0;

/////////////////////////////////////////////////////////////////////////////
//
// Structs & typedefs
//
/////////////////////////////////////////////////////////////////////////////

#define TA1 TNT::Array1D
#define TA2 TNT::Array2D

// Typical XYZ point and array of points
typedef struct {
    float x;
    float y;
    float z;
} cePoint, *pcePoint;

// An AFP (aligned fragment pair), and list/pointer
typedef struct {
    int first;
    int second;
} afp, *path, **pathCache;


// calculates a simple distance matrix
float** calcDM(const vector<vector<float> > &coords, int len)
{
    int i = 0;

    float** dm = (float**) malloc(sizeof(float*)*len);
    for ( i = 0; i < len; i++ )
        dm[i] = (float*) malloc( sizeof(float)*len);

    int row=0, col=0;
    for ( row = 0; row < len; row++ ) 
    {
        for ( col = 0; col < len; col++ ) 
        {
            dm[row][col]=sqrt(pow(coords[row][0] - coords[col][0],2) +
                              pow(coords[row][1] - coords[col][1],2) +
                              pow(coords[row][2] - coords[col][2],2) );
        }
    }
    return dm;
}

// Calculates the CE Similarity Matrix
float** calcS(float** d1, float** d2, int lenA, int lenB, float winSize)
{
    int i;
    
    // initialize the 2D similarity matrix
    float** S = (float**) malloc(sizeof(float*)*lenA);
    for ( i = 0; i < lenA; i++ )
        S[i] = (float*) malloc( sizeof(float)*lenB);

    float sumSize = (winSize-1.0)*(winSize-2.0) / 2.0;
    //
    // This is where the magic of CE comes out.  In the similarity matrix,
    // for each i and j, the value of ceSIM[i][j] is how well the residues
    // i - i+winSize in protein A, match to residues j - j+winSize in protein
    // B.  A value of 0 means absolute match; a value >> 1 means bad match.
    //
    int iA, iB, row, col;
    for ( iA = 0; iA < lenA; iA++ ) {
        for ( iB = 0; iB < lenB; iB++ ) {
            S[iA][iB] = -1.0;
            if ( iA > lenA - winSize || iB > lenB - winSize )
                continue;
        
            float score = 0.0;

            //
            // We always skip the calculation of the distance from THIS
            // residues, to the next residue.  This is a time-saving heur-
            // istic decision.  Almost all alpha carbon bonds of neighboring
            // residues is 3.8 Angstroms.  Due to entropy, S = -k ln pi * pi,
            // this tell us nothing, so it doesn't help so ignore it.
            //
            for ( row = 0; row < (int) winSize - 2; row++ ) {
                for ( col = row + 2; col < (int) winSize; col++ ) {
                    score += fabs( d1[iA+row][iA+col] - d2[iB+row][iB+col] );
                }
            }

            S[iA][iB] = score / sumSize;
        }
    }
    return S;
}



/* Converter: ChainUnit -> vector<vector<float> > */
vector<vector<float> > getCoords(const ChainUnit &chain)
{
    // make space for the current coords
    vector<float> tmp_array3(3,0.);
    int L=chain.residues.size();
    vector<vector<float> > coords(L,tmp_array3);

    // loop through the arguments, pulling out the XYZ coordinates.
    int r,a;
    for ( r=0; r<L; r++ )
    {
        for (a=0; a<chain.residues[r].atoms.size(); a++)
        {
            if (chain.residues[r].atoms[a].name==" CA ")
            {
                coords[r][0]=chain.residues[r].atoms[a].xyz[0];
                coords[r][1]=chain.residues[r].atoms[a].xyz[1];
                coords[r][2]=chain.residues[r].atoms[a].xyz[2];
            }
        }
    }
    return coords;
}


/* Optimal path finding algorithm (CE). */
pathCache findPath( float** S, float** dA, float** dB, int lenA, int lenB, int winSize, int& bufferSize )
{
    // CE-specific cutoffs
    const float D0 = 3.0;
    const float D1 = 4.0;
    const int MAX_KEPT = 20;
    const int gapMax = 30;

    // the best Path's score
    float bestPathScore = 1e6;
    int bestPathLength = 0;

    // length of longest possible alignment
    int smaller = ( lenA < lenB ) ? lenA : lenB;
    int winSum = (winSize-1)*(winSize-2)/2;

    //
    // BEST PATH
    //
    path bestPath = (path) malloc( sizeof(afp)*smaller );

    // index variable for below
    int i, j;
    for ( i = 0; i < smaller; i++ ) {
        bestPath[i].first = -1;
        bestPath[i].second = -1;
    }

    //======================================================================
    // for storing the best 20 paths
    int bufferIndex = 0; //, bufferSize = 0;
    int lenBuffer[MAX_KEPT];
    float scoreBuffer[MAX_KEPT];
    pathCache pathBuffer = (pathCache) malloc(sizeof(path*)*MAX_KEPT);

    for ( i = 0; i < MAX_KEPT; i++ ) {
        // initialize the paths
        scoreBuffer[i] = 1e6;
        lenBuffer[i] = 0;
        pathBuffer[i] = 0;
    }

    // winCache
    // this array stores a list of residues seen.  We use it to calculate the
    // total score of a path from 1..M and then add it to M+1..N.
    int* winCache = (int*) malloc(sizeof(int)*smaller);
    for ( i = 0; i < smaller; i++ )
        winCache[i] = (i+1)*i*winSize/2 + (i+1)*winSum;

    // allScoreBuffer
    // this 2D array keeps track of all partial gapped scores
    float** allScoreBuffer = (float**) malloc(sizeof(float*)*smaller);
    for ( i = 0; i < smaller; i++ ) {
        allScoreBuffer[i] = (float*) malloc( (gapMax*2+1) * sizeof(float));
        // initialize the ASB
        for ( j = 0; j < gapMax*2+1; j++ )
            allScoreBuffer[i][j] = 1e6;
    }
    
    int* tIndex = (int*) malloc(sizeof(int)*smaller);
    int gapBestIndex = -1;

    //======================================================================
    // Start the search through the CE matrix.
    //
    int iA, iB;
    for ( iA = 0; iA < lenA; iA++ ) {
        if ( iA > lenA - winSize*(bestPathLength-1) )
            break;

        for ( iB = 0; iB < lenB; iB++ ) {
            if ( S[iA][iB] >= D0 )
                continue;
            
            if ( S[iA][iB] == -1.0 )
                continue;
            
            if ( iB > lenB - winSize*(bestPathLength-1) )
                break;

            //
            // Restart curPath here.
            //
            path curPath = (path) malloc( sizeof(afp)*smaller );

            int i;
            for ( i = 0; i < smaller; i++ ) {
                curPath[i].first = -1;
                curPath[i].second = -1;
            }
            curPath[0].first = iA;
            curPath[0].second = iB;
            int curPathLength = 1;
            tIndex[curPathLength-1] = 0;
            float curTotalScore = 0.0;

            //
            // Check all possible paths starting from iA, iB
            //
            int done = 0;
            while ( ! done ) {
                float gapBestScore = 1e6;
                gapBestIndex = -1;
                int g;

                //
                // Check all possible gaps [1..gapMax] from here
                //
                for ( g = 0; g < (gapMax*2)+1; g++ ) {
                    int jA = curPath[curPathLength-1].first + winSize;
                    int jB = curPath[curPathLength-1].second + winSize;

                    if ( (g+1) % 2 == 0 ) {
                        jA += (g+1)/2;
                    }
                    else { // ( g odd )
                        jB += (g+1)/2;
                    }

                    //
                    // Following are three heuristics to ensure high quality
                    // long paths and make sure we don't run over the end of
                    // the S, matrix.
                    
                    // 1st: If jA and jB are at the end of the matrix
                    if ( jA > lenA-winSize-1 || jB > lenB-winSize-1 ){
                        continue;
                    }
                    // 2nd: If this gapped octapeptide is bad, ignore it.
                    if ( S[jA][jB] > D0 )
                        continue;
                    // 3rd: if too close to end, ignore it.
                    if ( S[jA][jB] == -1.0 )
                        continue;
                    
                    float curScore = 0.0;
                    int s;
                    for ( s = 0; s < curPathLength; s++ ) {
                        curScore += fabs( dA[curPath[s].first][jA] - dB[curPath[s].second][jB] );
                        curScore += fabs( dA[curPath[s].first  + (winSize-1)][jA+(winSize-1)] - 
                                          dB[curPath[s].second + (winSize-1)][jB+(winSize-1)] );
                        int k;
                        for ( k = 1; k < winSize-1; k++ )
                            curScore += fabs( dA[curPath[s].first  + k][ jA + (winSize-1) - k ] - 
                                              dB[curPath[s].second + k][ jB + (winSize-1) - k ] );
                    }
                    
                    curScore /= (float) winSize * (float) curPathLength;

                    if ( curScore >= D1 ) {
                        continue;
                    }

                    // store GAPPED best                    
                    if ( curScore < gapBestScore ) {
                        curPath[curPathLength].first = jA;
                        curPath[curPathLength].second = jB;
                        gapBestScore = curScore;
                        gapBestIndex = g;
                        allScoreBuffer[curPathLength-1][g] = curScore;
                    }
                } /// ROF -- END GAP SEARCHING
                
                    //
                    // DONE GAPPING:
                    //

                    // calculate curTotalScore
                    curTotalScore = 0.0;
                    int jGap, gA, gB;
                    float score1=0.0, score2=0.0;
                
                    if ( gapBestIndex != -1 ) {
                        jGap = (gapBestIndex + 1 ) / 2;
                        if ((gapBestIndex + 1 ) % 2 == 0) {
                            gA = curPath[ curPathLength-1 ].first + winSize + jGap;
                            gB = curPath[ curPathLength-1 ].second + winSize;
                        }
                        else {
                            gA = curPath[ curPathLength-1 ].first + winSize;
                            gB = curPath[ curPathLength-1 ].second + winSize + jGap;
                        }

                        // perfect
                        score1 = (allScoreBuffer[curPathLength-1][gapBestIndex] * winSize * curPathLength
                                  + S[gA][gB]*winSum)/(winSize*curPathLength+winSum);

                        // perfect
                        score2 = ((curPathLength > 1 ? (allScoreBuffer[curPathLength-2][tIndex[curPathLength-1]])
                                   : S[iA][iB])
                                  * winCache[curPathLength-1] 
                                  + score1 * (winCache[curPathLength] - winCache[curPathLength-1]))
                            / winCache[curPathLength];

                        curTotalScore = score2;
                        // heuristic -- path is getting sloppy, stop looking
                        if ( curTotalScore > D1 ) {
                            done = 1;
                            gapBestIndex=-1;
                            break;
                        }
                        else {
                            allScoreBuffer[curPathLength-1][gapBestIndex] = curTotalScore;
                            tIndex[curPathLength] = gapBestIndex;
                            curPathLength++;
                        }
                    }
                    else {
                        // if here, then there was no good gapped path
                        // so quit and restart from iA, iB+1
                        done = 1;
                        curPathLength--;
                        break;
                    }

                    //
                    // test this gapped path against the best seen
                    // starting from iA, iB
                    // 

                    // if our currently best gapped path from iA and iB is LONGER
                    // than the current best; or, it's equal length and the score's
                    // better, keep the new path.
                    if ( curPathLength > bestPathLength ||
                         (curPathLength == bestPathLength && curTotalScore < bestPathScore )) {
                        bestPathLength = curPathLength;
                        bestPathScore = curTotalScore;
                        // deep copy curPath
                        path tempPath = (path) malloc( sizeof(afp)*smaller );

                        int i;
                        for ( i = 0; i < smaller; i++ ) {
                            tempPath[i].first = curPath[i].first;
                            tempPath[i].second = curPath[i].second;
                        }

                        free(bestPath);
                        bestPath = tempPath;
                    }
            } /// END WHILE

                //
                // At this point, we've found the best path starting at iA, iB.
                //
                if ( bestPathLength > lenBuffer[bufferIndex] ||
                     ( bestPathLength == lenBuffer[bufferIndex] &&
                       bestPathScore < scoreBuffer[bufferIndex] )) {

                    // we're going to add an entry to the ring-buffer.
                    // Adjust maxSize values and curIndex accordingly.
                    bufferIndex = ( bufferIndex == MAX_KEPT-1 ) ? 0 : bufferIndex+1;
                    bufferSize = ( bufferSize < MAX_KEPT ) ? bufferSize+1 : MAX_KEPT;
                    path pathCopy = (path) malloc( sizeof(afp)*smaller );

                    int i;
                    for ( i = 0; i < smaller; i++ ) {
                        pathCopy[i].first = bestPath[i].first;
                        pathCopy[i].second = bestPath[i].second;
                    }

                    if ( bufferIndex == 0 && bufferSize == MAX_KEPT ) {
                        if ( pathBuffer[MAX_KEPT-1] )
                            free(pathBuffer[MAX_KEPT-1]); 
                        pathBuffer[MAX_KEPT-1] = pathCopy;
                        scoreBuffer[MAX_KEPT-1] = bestPathScore;
                        lenBuffer[MAX_KEPT-1] = bestPathLength;
                    }
                    else {    
                        if ( pathBuffer[bufferIndex-1] )
                            free(pathBuffer[bufferIndex-1]);
                        pathBuffer[bufferIndex-1] = pathCopy;
                        scoreBuffer[bufferIndex-1] = bestPathScore;
                        lenBuffer[bufferIndex-1] = bestPathLength;
                    }
                }
                free(curPath); curPath=0;
        } // ROF -- end for iB
    } // ROF -- end for iA

    return pathBuffer;
}

/* tranpose of a matrix */
TA2<float> transpose( const TA2<float>& v )
{
    unsigned int m = v.dim1();
    unsigned int n = v.dim2();
    
    TA2<float> rVal(n,m);
    
    for ( unsigned int i = 0; i < m; i++ )
        for ( unsigned int j = 0; j < n; j++ )
            rVal[j][i] = v[i][j];
    return rVal;
}

/* filter through the results and find the best */
void findBest(vector<vector<float> > &rVal, 
    const vector<vector<float> >& coordsA,
    const vector<vector<float> >& coordsB,
    pathCache paths, int bufferSize, int smaller, int winSize )
{
    // keep the best values
    float bestRMSD = 1e6;
    TA2<float> bestU;
    TA1<float> bestCOM1, bestCOM2;
    int bestLen = 0;

    // loop through the buffer
    for ( int o = 0; o < bufferSize; o++ ) 
    {
        // grab the current path
        TA2<float> c1(smaller, 3, 0.0);
        TA2<float> c2(smaller, 3, 0.0);
        int curLen = 0;

        int j = 0; int it = 0;
        while ( j < smaller ) 
        {
            // rebuild the coordinate lists for this path
            if ( paths[o][j].first != -1 )
            {
                int k = 0;
                while ( k++ < winSize )
                {
                    float t1[] = { coordsA[paths[o][j].first +k][0], 
                                    coordsA[paths[o][j].first +k][1], 
                                    coordsA[paths[o][j].first +k][2]};
                    float t2[] = { coordsB[paths[o][j].second+k][0],
                                    coordsB[paths[o][j].second+k][1],
                                    coordsB[paths[o][j].second+k][2]};
                    for ( int d = 0; d < c1.dim2(); d++ ) 
                    {
                        c1[it][d] =  t1[d];
                        c2[it][d] =  t2[d];
                    }
                    it++;
                }
                j++;
            }
            else 
            {
                curLen = it;
                break;
            }
        }

        // For convenience, let there be M points of N dimensions
        unsigned int m = curLen;
        unsigned int n = c2.dim2();

        /* Superpose the two proteins */

        // centers of mass for c1 and c2
        TA1<float> c1COM(n,0.0);
        TA1<float> c2COM(n,0.0); 
       
        // Calc CsOM
        for (unsigned int i = 0; i < m; i++ )
        {
            for (unsigned int j = 0; j < n; j++ )
            {
                c1COM[j] += (float) c1[i][j] / (float) m;
                c2COM[j] += (float) c2[i][j] / (float) m;
            }
        }

        // Move the two vectors to the origin
        for (unsigned int i = 0; i < m; i++ )
        {
            for (unsigned int j = 0; j < n; j++ )
            {
                c1[i][j] -= c1COM[j];
                c2[i][j] -= c2COM[j];
            }
        }

        //==========================================================================
        //    
        // Calculate U and RMSD.  This is broken down to the super-silly-easy    
        // math of: U = Wt * V, where Wt and V are NxN matrices from the SVD of     
        // R, the correlation matrix between the two origin-based vector sets.    
        //    
        //==========================================================================    
        
        // Calculate the initial residual, E0    
        // E0 = sum( Yn*Yn + Xn*Xn ) -- sum of squares    
        float E0 = 0.0;
        for (unsigned int i = 0; i < m; i++ )
        {
            for (unsigned int j = 0; j < n; j++ )
            {
                E0 += (c1[i][j]*c1[i][j])+(c2[i][j]*c2[i][j]);
            }
        }
        
        //        
        // SVD is the SVD of the correlation matrix Xt*Y    
        // R = c2' * c1 = W * S * Vt    
        JAMA::SVD<float> svd = JAMA::SVD<float>( TNT::matmult(transpose(c2), c1 ) );
    
        // left singular vectors
        TA2<float> W = TA2<float>(n,n);
        // diago    nal matrix of singular values
        TA2<float> S = TA2<float>(n,n);
        // right singular ve    ctors
        TA2<float> Vt = TA2<float>(n,n);
        // singular values        
        TA1<float> sigmas = TA1<float>(n);
        
        svd.getU(W);
        svd.getS(S);
        svd.getV(Vt);
        Vt = transpose(Vt);
        svd.getSingularValues(sigmas);
        
        //    
        // Check any reflections before rotation of the points;            
        // if det(W)*det(V) == -1 then we just reflect        
        // the principal axis corresponding to the smallest eigenvalue by -1    
        //        
        JAMA::LU<float> LU_Vt(Vt);
        JAMA::LU<float> LU_W(W);
        
        if ( LU_W.det() * LU_Vt.det() == -1 )
        {
            //std::cout << "_________REFLECTION_________" << std::endl;    
            
            // revese the smallest axes and last sigma    
            S[n-1][n-1] = -S[n-1][n-1];
            
            for ( unsigned int i = 0; i < n; i++ )
                W[n-1][i] = -W[n-1][i];
            
            sigmas[n-1] = -sigmas[n-1];
        }
        
        // calculate the rotation matrix, U.    
        // U = W * Vt    
        TA2<float> U = TA2<float>(TNT::matmult(W, Vt));

        //    
        // Rotate the points in the 2nd vector by U, thus solving the problem    
        //    
        c2 = TNT::matmult( c2, U );

        //    
        // Now calculate the RMSD    
        //    
        float sig = 0.0;
        for ( unsigned int i = 0; i < n; i++ )
            sig += sigmas[i];
        
        float curRMSD = sqrt(fabs((E0 - 2*sig) / (float) m ));

        //
        // Save the best
        //
        if ( curRMSD < bestRMSD || ( curRMSD == bestRMSD && c1.dim1() > bestLen )) {
            bestU = U.copy();
            bestRMSD = curRMSD;
            bestCOM1 = c1COM.copy();
            bestCOM2 = c2COM.copy();
            bestLen = curLen;
        }

        o++;
    }

    if ( bestRMSD == 1e6 ) {
        cerr << "ERROR: Best RMSD found was 1e6.  Broken.\n";
        return;
    }
    
    // list of list of pairs    
    rVal[0][0]=bestU[0][0]; rVal[0][1]=bestU[1][0]; rVal[0][2]=bestU[2][0];
    rVal[1][0]=bestU[0][1]; rVal[1][1]=bestU[1][1]; rVal[1][2]=bestU[2][1];
    rVal[2][0]=bestU[0][2]; rVal[2][1]=bestU[1][2]; rVal[2][2]=bestU[2][2];

    rVal[0][3]=bestCOM1[0];
    rVal[1][3]=bestCOM1[1];
    rVal[2][3]=bestCOM1[2];

    rVal[3][0]=-bestCOM2[0];
    rVal[3][1]=-bestCOM2[1];
    rVal[3][2]=-bestCOM2[2];

    rVal[4][0]=bestLen;
    rVal[4][1]=bestRMSD;
    return;
}

/* C-Cealign; driver. */
vector<vector<float> > ccealign_ccealign(
    const ChainUnit &chain_moving, const ChainUnit &chain_fixed)
{
    vector<float> tmp_array4(4,0.);
    vector<vector<float> >rVal(5,tmp_array4);

    // index variables
    int i=0;
    
    // handle empty selections (should probably do this in Python)
    int lenA = chain_fixed.residues.size();
    int lenB = chain_moving.residues.size();
    if ( lenB < 1 ) 
    {
        cerr<<"CEALIGN ERROR: Second selection didn't have any atoms.  Please check your selection.\n"<<endl;
        cerr<<"CEALIGN ERROR: Second selection didn't have any atoms.  Please check your selection.\n"<<endl;
        return rVal;
    }
    int smaller = lenA < lenB ? lenA : lenB;
    
    // get the coodinates from the Python objects
    vector<vector<float> > coordsA=getCoords(chain_fixed);
    vector<vector<float> > coordsB=getCoords(chain_moving);

    // calculate the distance matrix for each protein
    float** dmA = calcDM(coordsA, lenA);
    float** dmB = calcDM(coordsB, lenB);

    // calculate the CE Similarity matrix
    float **S = calcS(dmA, dmB, lenA, lenB, windowSize);

    // find the best path through the CE Sim. matrix
    int bufferSize = 0;
    int& bInt = bufferSize;
    pathCache paths = findPath( S, dmA, dmB, lenA, lenB, (int) windowSize, bInt );

    // Get the optimal superposition here...
    findBest(rVal, coordsA, coordsB, paths, bufferSize, smaller, (int) windowSize );

    // release memory
    coordsA.clear();
    coordsB.clear();

    // distance matrices
    for ( i = 0; i < lenA; i++ ) free( dmA[i] );
    free(dmA);

    for  ( i = 0; i < lenB; i++ ) free( dmB[i] );
    free(dmB);

    // similarity matrix
    for ( i = 0; i < lenA; i++ ) free( S[i] );
    free(S);

    return rVal;
}

#endif
