#include "PDBParser.hpp"

#include "basic_define.h"
#include "global_var.h"
#include "param_set.h"
#include "basic_fun.h"
#include "NW.h"
#include "Kabsch.h"
#include "TMalign.h"

using namespace std;

/* convert ChainUnit to array */
int read_PDB(const ChainUnit & chain, double **coor, char *seq, 
    int *resno, int **nres)
{
    int i=0;
    string line, str;    
    string du1, i8;
    
    int mk = 1;
    int r,a;
    for (r=0;r<chain.residues.size();r++)
    {
        if (chain.residues[r].het==true) continue;
        du1 = chain.residues[r].icode;
        int nDu1 = *(du1.c_str());// the ASCII code of du1

        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name!=" CA ") continue;

            //get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]); 
            for (int j=0;j<3;j++) coor[i][j]=chain.residues[r].atoms[a].xyz[j];
            seq[i]=aa3to1(chain.residues[r].resn);
            resno[i]=chain.residues[r].resi;

            nres[resno[i]][nDu1] ++;
            i++;
        }
    }
    seq[i]='\0';   
    
    return i;
}

/* allocate memory according to input ChainUnit */
void load_PDB_allocate_memory(const ChainUnit & chain1, 
    const ChainUnit & chain2)
{    
    /* load data */
    NewArray(&nres1, NMAX, ASCIILimit);// Only data from nres1[0~ NMAX-1][32~122] is used
    NewArray(&nres2, NMAX, ASCIILimit);
    for (int i = 0; i < NMAX; i++)// Initialization
    {
        for (int j = 0; j < ASCIILimit; j++)
        {
            nres1[i][j] = 0;
            nres2[i][j] = 0;
        }
    }

    tempxlen = chain1.residues.size();// Get predicted length
    tempylen = chain2.residues.size();

    /* allocate memory for x and y */
    NewArray(&xa, tempxlen, 3);
    seqx = new char[tempxlen + 1];
    secx = new int[tempxlen];
    xresno = new int[tempxlen]; 
   
    NewArray(&ya, tempylen, 3);
    seqy = new char[tempylen + 1];
    yresno = new int[tempylen];
    secy = new int[tempylen];

    xlen = read_PDB(chain1, xa, seqx, xresno, nres1);// Get exact length
    ylen = read_PDB(chain2, ya, seqy, yresno, nres2);
    minlen = min(xlen, ylen);
    
    /* allocate memory for other temporary varialbes */
    NewArray(&r1, minlen, 3);
    NewArray(&r2, minlen, 3);
    NewArray(&xtm, minlen, 3);
    NewArray(&ytm, minlen, 3);
    NewArray(&xt, xlen, 3);

    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);  
}

/* stick to alignment. for some reason, its TM-score is less than the 
 * TM-score of Jianyi's TMalignc */
void TMalign_I(const string & aln1, const string & aln2,
    ChainUnit & chain1, const ChainUnit & chain2,
    float & tmscore1, float & tmscore2, int f=0)
{
    tmscore1=tmscore2=0;

    fast_level = 0;
    f_opt = false;
    if (f>0)
    {
        fast_level=f;
        f_opt=true;
    }

    /* read initial alignment */
    strcpy(sequence[0],aln1.c_str());
    strcpy(sequence[1],aln2.c_str());

    /* convert data */ 
    load_PDB_allocate_memory(chain1, chain2);

    /* parameter set */
    parameter_set4search(xlen, ylen);//please set parameters in the function
    int simplify_step     = 40;      //for similified search engine
    int score_sum_method  = 8;       //for scoring method, whether only sum
                                     //over pairs with dis<score_d8
        
    int i;
    int *invmap0          = new int[ylen+1]; 
    int *invmap           = new int[ylen+1]; 
    double TM, TMmax=-1;
    for(i=0; i<ylen; i++) invmap0[i]=-1;


    double ddcc=0.4;
    if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    double local_d0_search = d0_search;

    /* parse initial alignment */
    for (int j = 1; j < ylen; j++)// Set aligned position to be "-1"
        invmap[j] = -1;

    int i1 = -1;// in C version, index starts from zero, not from one
    int i2 = -1;
    int L1 = strlen(sequence[0]);
    int L2 = strlen(sequence[1]);
    int L = min(L1, L2);// Get positions for aligned residues
    for (int kk1 = 0; kk1 < L; kk1++)
    {
        if (sequence[0][kk1] != '-') i1++;
        if (sequence[1][kk1] != '-')
        {
            i2++;
            if (i2 >= ylen || i1 >= xlen) kk1 = L;
            else if (sequence[0][kk1] != '-') invmap[i2] = i1;
        }
    }

    /* 2. Align proteins from original alignment */
    double prevD0_MIN = D0_MIN;// stored for later use
    int prevLnorm = Lnorm;
    double prevd0 = d0;
    TM_ali = standard_TMscore(xa, ya, xlen, ylen, invmap, L_ali, rmsd_ali);
    D0_MIN = prevD0_MIN;
    Lnorm = prevLnorm;
    d0 = prevd0;
    TM = detailed_search_standard(xa, ya, xlen, ylen, invmap, t, u, 40, 8, 
        local_d0_search, true);
    if (TM > TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }


    /* The alignment will not be changed any more in the following */
    bool flag=false;
    for(i=0; i<ylen; i++)
    {
        if(invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag) return;


    /* Detailed TMscore search engine --> prepare for final TMscore */
    /* run detailed TMscore search engine for the best alignment, and */
    /* extract the best rotation matrix (t, u) for the best alginment */
    simplify_step=1;
    score_sum_method=8;
    TM = detailed_search_standard(xa, ya, xlen, ylen, invmap0, t, u,
        simplify_step, score_sum_method, local_d0_search, false);


    /* select pairs with dis<d8 for final TMscore computation and output alignment */
    int n_ali8, k=0;
    int n_ali=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    do_rotation(xa, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;        
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
            if (d <= score_d8 || (I_opt == true))
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];
                    
                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2]; 
               
                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya[j][0];
                r2[k][1] = ya[j][1];
                r2[k][2] = ya[j][2]; 
               
                k++;
            }
        }
    }
    n_ali8=k;

    double rmsd0 = 0.0;
    Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);//rmsd0 is used for final output, 
                                            //only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);

    /* Final TMscore */
    double rmsd;
    double t0[3], u0[3][3];
    double TM1, TM2;
    double d0_out=5.0;  
    simplify_step=1;
    score_sum_method=0;

    double d0_0, TM_0;
    double Lnorm_0=ylen;


    /* normalized by length of structure A */
    parameter_set4final(Lnorm_0);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
    TM_0 = TM1;

    /* normalized by length of structure B */
    parameter_set4final(xlen+0.0);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);

    /* Done! Free memory */
    free_memory();

    delete [] invmap0;
    delete [] invmap;
    delete [] m1;
    delete [] m2;
    sequence[0][0]=0;
    sequence[0][1]=0;

    /* return value */
    tmscore2=TM1;
    tmscore1=TM2;
    return;
}


/* free alignment */
void TMalign(string & aln1, string & aln2,
    ChainUnit & chain1, const ChainUnit & chain2,
    float & tmscore1, float & tmscore2, int f=0)
{
    tmscore1=tmscore2=0;
    aln1.clear();
    aln2.clear();

    fast_level = 0;
    f_opt = false;
    if (f>0)
    {
        fast_level=f;
        f_opt=true;
    }

    /* convert data */ 
    load_PDB_allocate_memory(chain1, chain2);

    /* parameter set */
    parameter_set4search(xlen, ylen);//please set parameters in the function
    int simplify_step     = 40;      //for similified search engine
    int score_sum_method  = 8;       //for scoring method, whether only sum
                                     //over pairs with dis<score_d8
        
    int i;
    int *invmap0          = new int[ylen+1]; 
    int *invmap           = new int[ylen+1]; 
    double TM, TMmax=-1;
    for(i=0; i<ylen; i++) invmap0[i]=-1;


    double ddcc=0.4;
    if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
    double local_d0_search = d0_search;

    /* get initial alignment with gapless threading */
    get_initial(xa, ya, xlen, ylen, invmap0);
    //find the max TMscore for this initial alignment with the simplified search_engin
    TM = detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, 
        score_sum_method, local_d0_search);
    if (TM>TMmax) TMmax = TM;
    //run dynamic programing iteratively to find the best alignment
    TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }


    /* get initial alignment based on secondary structure */
    get_initial_ss(xa, ya, xlen, ylen, invmap);
    TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, 
        score_sum_method, local_d0_search);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*0.2)
    {
        TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /*    get initial alignment based on local superposition    */
    //=initial5 in original TM-align
    if (get_initial5(xa, ya, xlen, ylen, invmap))
    {
        TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, 
            score_sum_method, local_d0_search);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
        if (TM > TMmax*ddcc)
        {
            TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 2, 
                local_d0_search);
            if (TM>TMmax)
            {
                TMmax = TM;
                for (int i = 0; i<ylen; i++) invmap0[i] = invmap[i];
            }
        }
    }
    else
    {
        cout << "\n\nWarning: initial alignment from local "
            "superposition fail!\n\n";
    }


    /* get initial alignment based on previous alignment+secondary structure */
    //=initial3 in original TM-align
    get_initial_ssplus(xa, ya, xlen, ylen, invmap0, invmap);
    TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, 
        score_sum_method, local_d0_search);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /*    get initial alignment based on fragment gapless threading    */
    //=initial4 in original TM-align
    get_initial_fgt(xa, ya, xlen, ylen, xresno, yresno, invmap);
    TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (i = 0; i<ylen; i++)
        {
            invmap0[i] = invmap[i];
        }
    }
    if (TM > TMmax*ddcc)
    {
        TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (i = 0; i<ylen; i++) invmap0[i] = invmap[i];
        }
    }


    /* The alignment will not be changed any more in the following */
    bool flag=false;
    for(i=0; i<ylen; i++)
    {
        if(invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag) return;


    /* Detailed TMscore search engine --> prepare for final TMscore */
    /* run detailed TMscore search engine for the best alignment, and */
    /* extract the best rotation matrix (t, u) for the best alginment */
    simplify_step=1;
    score_sum_method=8;
    TM = detailed_search_standard(xa, ya, xlen, ylen, invmap0, t, u,
        simplify_step, score_sum_method, local_d0_search, false);


    /* select pairs with dis<d8 for final TMscore computation and output alignment */
    int n_ali8, k=0;
    int n_ali=0;
    int *m1, *m2;
    double d;
    m1=new int[xlen]; //alignd index in x
    m2=new int[ylen]; //alignd index in y
    do_rotation(xa, xt, xlen, t, u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            n_ali++;        
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
            if (d <= score_d8 || (I_opt == true))
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];
                    
                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2]; 
               
                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya[j][0];
                r2[k][1] = ya[j][1];
                r2[k][2] = ya[j][2]; 
               
                k++;
            }
        }
    }
    n_ali8=k;

    double rmsd0 = 0.0;
    Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);//rmsd0 is used for final output, 
                                            //only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);

    /* Final TMscore */
    double rmsd;
    double t0[3], u0[3][3];
    double TM1, TM2;
    double d0_out=5.0;  
    simplify_step=1;
    score_sum_method=0;

    double d0_0, TM_0;
    double Lnorm_0=ylen;


    /* normalized by length of structure A */
    parameter_set4final(Lnorm_0);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
    TM_0 = TM1;

    /* normalized by length of structure B */
    parameter_set4final(xlen+0.0);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);

    /* retrieve alignment */
    int ali_len=xlen+ylen; //maximum length of alignment
    char *seqxA, *seqyA;
    seqxA=new char[ali_len];
    seqyA=new char[ali_len];
    int kk=0, i_old=0, j_old=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        kk++;  
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }
    aln1=(string)(seqxA);
    aln2=(string)(seqyA);
    delete [] seqxA;
    delete [] seqyA;

    /* Done! Free memory */
    free_memory();

    delete [] invmap0;
    delete [] invmap;
    delete [] m1;
    delete [] m2;

    /* return value */
    tmscore2=TM1;
    tmscore1=TM2;
    return;
}
