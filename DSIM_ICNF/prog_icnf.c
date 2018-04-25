
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

#define MAIN
#include "decl_ilib.h"
#include "decl_ivar.h"
#include "decl_ifun.h"

int main(int argc, char **argv)
{
    char comd[100];
    time_t iclk;
    FILE *fptr;

    if (argc != 6 && argc != 7)
    {
        fprintf(stdout,
            "This utility requires three (3) arguments:\n");
        fprintf(stdout, "[1] sequence input file\n");
        fprintf(stdout, "[2] DNA type [0 - B-DNA, 1 - Curved B-DNA, 2 - A-DNA]\n");
        fprintf(stdout, "[3] switch for complementarity\n");
        fprintf(stdout, "[4] directory name to store output\n");
        fprintf(stdout, "[5] ions flag\n");
        fprintf(stdout, "[6] if ions flag is not zero, ions input file\n");
        fprintf(stdout, "   if ions flag == 1 make DNA with ions\n");
        fprintf(stdout, "   if ions flag == 2 make DNA empty box\n");
        fprintf(stdout, "Please try again.\n");
        exit(1);
    }

    ions_flag = atoi(argv[5]);
    dna_type = atoi(argv[2]);

    // A bit of bulletproofing
    if (NULL == (fptr = fopen(argv[1], "r")))
    {
        fprintf(stdout, "Input file %s is not accessible.\n",
            argv[1]);
        exit(1);
    }
    fclose(fptr);

    if (ions_flag)
    {
        if (NULL == (fptr = fopen(argv[6], "r")))
        {
            fprintf(stdout, "Input file %s is not accessible.\n",
                argv[6]);
            exit(1);
        }
        fclose(fptr);
    }

    mkdir(argv[4], S_IRWXU);
    comp = atoi(argv[3]);
    time(&iclk);

    if ((dna_type == 0) || (dna_type == 1)) {
        genr_rtmp_bdna();
    } else if (dna_type == 2) {
         genr_rtmp_adna();
    } else {
        fprintf(stdout, "I don't recognize the DNA type %s\n",
            argv[2]);
        exit(1);
    }

    if (ions_flag == 2)
    {
        site.dna_totl = 0;
        sequ.nste = 0;
        seqv.nste = 0;
    } else {
        read_sequ(argv[1]);
    }
    if (ions_flag)
    {
        read_ions(argv[6]);
        sprintf(comd, "cp -f %s %s", argv[6], argv[4]);
        system(comd);
    }
    sprintf(comd, "cp -f %s %s", argv[1], argv[4]);
    system(comd);

    genr_dnac();
    if (ions_flag == 2)
    {
        cbon = 0;
        cben = 0;
        ctor = 0;
    } else {
        genr_topo();
    }
    if (ions_flag)
    {
        genr_ions();
    }
    wrte_rpsf(argv[4]);
    wrte_rxyz(argv[4]);
    wrte_lammps(argv[4]);
    return (0);
}
