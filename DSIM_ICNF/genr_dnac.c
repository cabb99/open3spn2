
/*********************************************************************
**********************************************************************
****        <<<<< 3SPN.2 Configuration Generator    >>>>>         ****
****                                                              ****  
****        Dan Hinckley    `                                     ****
****        Instutite for Molecular Engineering                   ****
****        University of Chicago                                 ****
****        Chicago, IL 60637                                     ****
****                                                              ****
**********************************************************************
*********************************************************************/

#include "decl_ilib.h"
#include "decl_ivar.h"
#include "decl_ifun.h"

void genr_dnac(void)
{
    long i, j, k, icnt;
    double dtha, dzdf, tshf, zshf;

    j = site.dna2 - 1;

    dtha = twist * _PI_ / 180.0;
    dzdf = rise;
    icnt = 0;

    for (i = sequ.nste - 1; i >= 0; i--)
    {
        tshf = dtha * (double)icnt;
        zshf = dzdf * (double)icnt;

        switch (sequ.base[i])
        {
            case 'A':
                atom[j].xval = A.rval * cos(A.thta + tshf);
                atom[j].yval = A.rval * sin(A.thta + tshf);
                atom[j].zval = A.zval + zshf;
                atom[j].mass = A.mass;
                sprintf(atom[j].snme, "A");
                sprintf(atom[j].rnme, "ADE");
                sprintf(atom[j - 1].rnme, "ADE");
                if (i)
                {
                    sprintf(atom[j - 2].rnme, "ADE");
                }
                atom[j].rnum = i;
                atom[j].stid = A.stid;
                atom[j].stid2 = A.stid2;
                atom[j].chrg = A.chrg;
                j = j - 1;
                break;
            case 'G':
                atom[j].xval = G.rval * cos(G.thta + tshf);
                atom[j].yval = G.rval * sin(G.thta + tshf);
                atom[j].zval = G.zval + zshf;
                atom[j].mass = G.mass;
                sprintf(atom[j].snme, "G");
                sprintf(atom[j].rnme, "GUA");
                sprintf(atom[j - 1].rnme, "GUA");
                atom[j].rnum = i;
                if (i)
                {
                    sprintf(atom[j - 2].rnme, "GUA");
                }
                atom[j].stid = G.stid;
                atom[j].stid2 = G.stid2;
                atom[j].chrg = G.chrg;
                j = j - 1;
                break;
            case 'T':
                atom[j].xval = T.rval * cos(T.thta + tshf);
                atom[j].yval = T.rval * sin(T.thta + tshf);
                atom[j].zval = T.zval + zshf;
                atom[j].mass = T.mass;
                sprintf(atom[j].snme, "T");
                sprintf(atom[j].rnme, "THY");
                sprintf(atom[j - 1].rnme, "THY");
                atom[j].rnum = i;
                if (i)
                {
                    sprintf(atom[j - 2].rnme, "THY");
                }
                atom[j].stid = T.stid;
                atom[j].stid2 = T.stid2;
                atom[j].chrg = T.chrg;
                j = j - 1;
                break;
            case 'C':
                atom[j].xval = C.rval * cos(C.thta + tshf);
                atom[j].yval = C.rval * sin(C.thta + tshf);
                atom[j].zval = C.zval + zshf;
                atom[j].mass = C.mass;
                sprintf(atom[j].snme, "C");
                sprintf(atom[j].rnme, "CYT");
                sprintf(atom[j - 1].rnme, "CYT");
                atom[j].rnum = i;
                if (i)
                {
                    sprintf(atom[j - 2].rnme, "CYT");
                }
                atom[j].stid = C.stid;
                atom[j].stid2 = C.stid2;
                atom[j].chrg = C.chrg;
                j = j - 1;
                break;
        }

        atom[j].xval = S.rval * cos(S.thta + tshf);
        atom[j].yval = S.rval * sin(S.thta + tshf);
        atom[j].zval = S.zval + zshf;
        atom[j].mass = S.mass;
        sprintf(atom[j].snme, "S");
        atom[j].rnum = i;
        atom[j].stid = S.stid;
        atom[j].stid2 = S.stid;
        atom[j].chrg = S.chrg;
        j = j - 1;

        if (i)
        {
            atom[j].xval = P.rval * cos(P.thta + tshf);
            atom[j].yval = P.rval * sin(P.thta + tshf);
            atom[j].zval = P.zval + zshf;
            atom[j].mass = P.mass;
            sprintf(atom[j].snme, "P");
            atom[j].rnum = i;
            atom[j].stid = P.stid;
            atom[j].stid2 = P.stid;
            if (ions_flag)
            {
                atom[j].chrg = P.chrg2;
            }
            else
            {
                atom[j].chrg = P.chrg;
            }
            j = j - 1;
        }
        icnt = icnt + 1;
    }
    

    icnt = 0;
    if (comp)
    {
        k = site.dna2;
        for (i = sequ.nste - 1; i >= 0; i--)
        {

            tshf = dtha * (double)icnt;
            zshf = dzdf * (double)icnt;

            if (i != sequ.nste - 1)
            {
                atom[k].xval = P.rval * cos(-P.thta + tshf);
                atom[k].yval = P.rval * sin(-P.thta + tshf);
                atom[k].zval = -P.zval + zshf;
                atom[k].mass = P.mass;
                sprintf(atom[k].snme, "P");
                atom[k].rnum = icnt + sequ.nste;
                atom[k].stid = P.stid;
                atom[k].stid2 = P.stid;
                if (ions_flag)
                {
                    atom[k].chrg = P.chrg2;
                }
                else
                {
                    atom[k].chrg = P.chrg;
                }
                k = k + 1;
            }

            atom[k].xval = S.rval * cos(-S.thta + tshf);
            atom[k].yval = S.rval * sin(-S.thta + tshf);
            atom[k].zval = -S.zval + zshf;
            atom[k].mass = S.mass;
            sprintf(atom[k].snme, "S");
            atom[k].rnum = icnt + sequ.nste;
            atom[k].stid = S.stid;
            atom[k].stid2 = S.stid;
            atom[k].chrg = S.chrg;
            k = k + 1;

            if (seqv.nste == 0)
            {

                switch (sequ.base[i])
                {
                    case 'T':
                        atom[k].xval = A.rval * cos(-A.thta + tshf);
                        atom[k].yval = A.rval * sin(-A.thta + tshf);
                        atom[k].zval = -A.zval + zshf;
                        atom[k].mass = A.mass;
                        sprintf(atom[k].snme, "A");
                        sprintf(atom[k].rnme, "ADE");
                        sprintf(atom[k - 1].rnme, "ADE");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "ADE");
                        }
                        atom[k].stid = A.stid;
                        atom[k].stid2 = A.stid2;
                        atom[k].chrg = A.chrg;
                        k = k + 1;
                        break;
                    case 'C':
                        atom[k].xval = G.rval * cos(-G.thta + tshf);
                        atom[k].yval = G.rval * sin(-G.thta + tshf);
                        atom[k].zval = -G.zval + zshf;
                        atom[k].mass = G.mass;
                        sprintf(atom[k].snme, "G");
                        sprintf(atom[k].rnme, "GUA");
                        sprintf(atom[k - 1].rnme, "GUA");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "GUA");
                        }
                        atom[k].stid = G.stid;
                        atom[k].stid2 = G.stid2;
                        atom[k].chrg = G.chrg;
                        k = k + 1;
                        break;
                    case 'A':
                        atom[k].xval = T.rval * cos(-T.thta + tshf);
                        atom[k].yval = T.rval * sin(-T.thta + tshf);
                        atom[k].zval = -T.zval + zshf;
                        atom[k].mass = T.mass;
                        sprintf(atom[k].snme, "T");
                        sprintf(atom[k].rnme, "THY");
                        sprintf(atom[k - 1].rnme, "THY");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "THY");
                        }
                        atom[k].stid = T.stid;
                        atom[k].stid2 = T.stid2;
                        atom[k].chrg = T.chrg;
                        k = k + 1;
                        break;
                    case 'G':
                        atom[k].xval = C.rval * cos(-C.thta + tshf);
                        atom[k].yval = C.rval * sin(-C.thta + tshf);
                        atom[k].zval = -C.zval + zshf;
                        atom[k].mass = C.mass;
                        sprintf(atom[k].snme, "C");
                        sprintf(atom[k].rnme, "CYT");
                        sprintf(atom[k - 1].rnme, "CYT");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "CYT");
                        }
                        atom[k].stid = C.stid;
                        atom[k].stid2 = C.stid2;
                        atom[k].chrg = C.chrg;
                        k = k + 1;
                        break;
                }
            }
            else
            {
                switch (seqv.base[i])
                {
                    case 'A':
                        atom[k].xval = A.rval * cos(-A.thta + tshf);
                        atom[k].yval = A.rval * sin(-A.thta + tshf);
                        atom[k].zval = -A.zval + zshf;
                        atom[k].mass = A.mass;
                        sprintf(atom[k].snme, "A");
                        sprintf(atom[k].rnme, "ADE");
                        sprintf(atom[k - 1].rnme, "ADE");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "ADE");
                        }
                        atom[k].stid = A.stid;
                        atom[k].stid2 = A.stid2;
                        atom[k].chrg = A.chrg;
                        k = k + 1;
                        break;
                    case 'G':
                        atom[k].xval = G.rval * cos(-G.thta + tshf);
                        atom[k].yval = G.rval * sin(-G.thta + tshf);
                        atom[k].zval = -G.zval + zshf;
                        atom[k].mass = G.mass;
                        sprintf(atom[k].snme, "G");
                        sprintf(atom[k].rnme, "GUA");
                        sprintf(atom[k - 1].rnme, "GUA");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "GUA");
                        }
                        atom[k].stid = G.stid;
                        atom[k].stid2 = G.stid2;
                        atom[k].chrg = G.chrg;
                        k = k + 1;
                        break;
                    case 'T':
                        atom[k].xval = T.rval * cos(-T.thta + tshf);
                        atom[k].yval = T.rval * sin(-T.thta + tshf);
                        atom[k].zval = -T.zval + zshf;
                        atom[k].mass = T.mass;
                        sprintf(atom[k].snme, "T");
                        sprintf(atom[k].rnme, "THY");
                        sprintf(atom[k - 1].rnme, "THY");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "THY");
                        }
                        atom[k].stid = T.stid;
                        atom[k].stid2 = T.stid2;
                        atom[k].chrg = T.chrg;
                        k = k + 1;
                        break;
                    case 'C':
                        atom[k].xval = C.rval * cos(-C.thta + tshf);
                        atom[k].yval = C.rval * sin(-C.thta + tshf);
                        atom[k].zval = -C.zval + zshf;
                        atom[k].mass = C.mass;
                        sprintf(atom[k].snme, "C");
                        sprintf(atom[k].rnme, "CYT");
                        sprintf(atom[k - 1].rnme, "CYT");
                        atom[k].rnum = icnt + sequ.nste;
                        if (i != sequ.nste - 1)
                        {
                            sprintf(atom[k - 2].rnme, "CYT");
                        }
                        atom[k].stid = C.stid;
                        atom[k].stid2 = C.stid2;
                        atom[k].chrg = C.chrg;
                        k = k + 1;
                        break;
                }
            }
            icnt = icnt + 1;
        }
    }
    genr_trns();
}
