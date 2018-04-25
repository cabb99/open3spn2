
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

#define _NA_ 6.022045E23
#define _PI_ M_PI
#define _TW_ 150

#ifndef MAIN
extern
#endif
    struct moie_data
{
    double xval;
    double yval;
    double zval;
    double rval;
    double thta;
    double mass;
    double chrg;
    double chrg2;
    long stid;
    long stid2;
} P, S, A, G, T, C, Na, Mg, Cl;

#ifndef MAIN
extern
#endif
    struct sequ_data
{
    long nste;
    char *base;
} sequ, seqv;

#ifndef MAIN
extern
#endif
    struct atom_data
{
    double xval;
    double yval;
    double zval;
    double mass;
    double chrg;
    char snme[4];
    char rnme[4];
    long rnum;
    long stid;
    long stid2;
} *atom, *ions;

#ifndef MAIN
extern
#endif
    struct site_data
{
    long dna_totl;
    long totl;
    long nions;
    long extra_Na;
    long extra_Cl;
    long dna1;
    long dna2;
} site;

#ifndef MAIN
extern
#endif
    struct resi_data
{
    long totl;
    long dna1;
    long dna2;
} resi;

#ifndef MAIN
extern
#endif
long cbon;

#ifndef MAIN
extern
#endif
    struct bond_data
{
    long aste;
    long bste;
} *bond;

#ifndef MAIN
extern
#endif
long cben, cben2;

#ifndef MAIN
extern
#endif
    struct bend_data
{
    long aste;
    long bste;
    long cste;
} *bend;

#ifndef MAIN
extern
#endif
long ctor;

#ifndef MAIN
extern
#endif
    struct tors_data
{
    long aste;
    long bste;
    long cste;
    long dste;
} *tors;

#ifndef MAIN
extern
#endif
    struct boxs_data
{
    double lenx;
    double leny;
    double lenz;
} boxs;

#ifndef MAIN
extern
#endif
long comp, idum, nIons, ions_flag, dna_type;
#ifndef MAIN
extern
#endif
double conc_MgCl2, conc_NaCl, twist, rise, ionBoxSize;

