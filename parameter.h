#ifndef PARAMETER_H
#define PARAMETER_H


typedef struct {

    /* phyical parameters */
    double BoxSize;
    double OmegaM0;
    double OmegaX0;
    double Hubble0;

    /* system parameters*/
    int NumSocket;
    int NumCpuPerSocket;
    int NumMicPerSocket;
    int NumThreadPerSocket;
    unsigned long int MemoBytesPerSocket;

    /* runtime parameters */
    int NumBits; // for mesh
    int PartBits; // for part
    int NumGridPerSide;
    unsigned long int TotNumPart;
    double OpenAngle;
    double CutRadius;
    double SoftenLength;
    char WorkPath[240];
    char SnapName[100];
    int  NumFile;

    /* status parameters */
    double redshift;

} GlobalParam;

void load_parameters(GlobalParam* gp, const char ParamFileName[]);
void setup_parameters(GlobalParam* gp);

#endif /* PARAMETER_H */

