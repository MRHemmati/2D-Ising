
/******************************************************************************/
/*** Computational Physics Course                                           ***/
/*** Name :   Sina Mehboodi & Mohammad Hemati & Nikta JabbarZade            ***/
/*** Stu Num: 941128        & 951133          & 951124                      ***/
/*** Ver 1.9.8                                                              ***/
/*** Date: 13960423                                                         ***/
/*** OS: Linux 16.10 (x86_64)                                               ***/
/*** Run under a Intel® Core™ i5-2450M CPU @ 2.50GHz - up to 3.1GHz         ***/
/*** 4GB DDR3 Memory                                                        ***/
/******************************************************************************/

// You can compile this code by following command in shell:
// > g++ -o Ising.out Ising.cpp -Ofast
// , then execute it in shell by
// > ./Ising.out

/* This is a sample code. It is not complete and it may contain some errors. 
 * So you can try it with your own risk! */

/* Random number generator in the following code is based on the random number generator function

   double ran2(long *idum);
   
implemented in the Numerical recipes in C, chapter 7 (ran2)

Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.

Visit www.nr.com for the licence.*/

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

using namespace std;

//BEGIN_FOLD - Random number generator section
//------------------------------------------------------------------------------

/* following routine based on Random number generator

   double ran2(long *idum);
   
implemented in the Numerical recipes in C, chapter 7 (ran2)

Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.

Visit www.nr.com for the licence.*/

// This is a internal, 32 bit random number generator with uniform Distribution in range [0..1)
/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum) {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

//------------------------------------------------------------------------------

// My interface 

// Random number generator seed
long iseed = -36;

// You can initializing random number generator seed by usin following routine
void Randomize() {
  iseed = -time(NULL);  
}

// return random number with uniform distribution in range [0..1).
inline double Random() { return ran2(&iseed); }

// return integer random number with uniform distribution in range 0 to N.
inline int Random(int N) { return int(ran2(&iseed)*N); }
//END_FOLD

const int L = 256;                // Lattice size; L = 2^n
const int LMask = L - 1;          // Mask used to perfom periodic boundary condition in a Perid() function
int** s;                          // Spins lattice; s[i][j] shows spin of cellm (i,j)
// All the following parameters are in the reduced unit system
float h = 0.0;                    // external field
float J = +1;                     // spin-spin coupling constant
float T = 2.2;                    // Tempreture
ofstream snapshot("snapshot.txt", // Keeps model snapshot
                  ios::out | ios::trunc);  
////////////////////////////////////////////////////////////////////////////////////////
ofstream MbeforT("MbeforT.txt");
ofstream magnTime15("magnTime15.txt");
ofstream magnTime20("magnTime20.txt");
ofstream magnTime225("magnTime225.txt");
ofstream magnTime40("magnTime40.txt");
ofstream magnTemp ("magnTemp.txt");
ofstream Susceptibility("Susceptibility.txt");
ofstream Corelationdata("Corelationdata.txt");
ofstream VarMagnetismPerT("VarMagnetismPerT.txt");
////////////////////////////////////////////////////////////////////////////////////////
double totalMmc=0;
double meanMmc=0;
int MCStep=0;
double meanMmcPerT=0;
double meanMmc2=0;
double totalMmc2=0;
double meanMmcPerT2=0;
double varM=0;
double var2M=0;
double Sus=0;

////////////////////////////////////////////////////////////////////////////////////////
void Initialize_text_file();
void Plot();
void Init();                      // Initialize model, alocate dynamic memory, ...
void Execute();                   // Simulate Ising model
void PostProcess();               // Perform post statistical analysis
void Done();                      // Dealocate dynamic memory, ...
/////////////////////////////////////////////////////////////////////////////////////////
//Plot function
void Plot(){
    cout<<"Plotting ....."<<endl;
    magnTemp.close();
    Susceptibility.close();
    magnTime15.close();
    magnTime20.close();
    magnTime225.close();
    magnTime40.close();
    Corelationdata.close();
    VarMagnetismPerT.close();
        cout<<"magnTemp"<<endl;
    system("gnuplot 'magnTemp.txt'");
        cout<<"energyTemp"<<endl;
        cout<<"Susceptibility"<<endl;
    system("gnuplot 'Susceptibility.txt'");
        cout<<"magnTime15"<<endl;
    system("gnuplot 'magnTime15.txt'");
        cout<<"magnTime20"<<endl;
    system("gnuplot 'magnTime20.txt'");
        cout<<"magnTime225"<<endl;
    system("gnuplot 'magnTime225.txt'");
        cout<<"magnTime40"<<endl;
    system("gnuplot 'magnTime40.txt'");
        cout<<"corelation"<<endl;
    system("gnuplot 'Corelationdata.txt'");
        cout<<"VarMagnetismPerT"<<endl;
    system("gnuplot 'VarMagnetismPerT.txt'");
}

//OUTput file
void Initialize_text_file(){

                   

           magnTemp<<"set title\"Variation of magnetism with respect to temperature\"\n"
                   <<"set terminal png \n"
                   <<"set output \"MagnPerTemp.png\"\n"
                   <<"set ylabel \"Magnetism\"\n"
                   <<"set xlabel\"Temperature\"\n"
                   <<"plot [][-.5:1.2] '-' using 1:2 title \"variation of magnetism\" w l lc 1"<<endl;
                
                           
Susceptibility<<"set title\"Variation of susceptibility with respect to temperature\"\n"
                           <<"set terminal png \n"
                           <<"set output \"Susceptibility.png\"\n"
                           <<"set ylabel \"Susceptibility\"\n"
                           <<"set xlabel\"Temperature\"\n"
                           <<"plot [][] '-' using 1:2 title \"variation of susceptibility\" w lp ps 1 pt 6"<<endl;
                           
magnTime15<<"set title\"Variation of magnetism with respect to time\"\n"
                   <<"set terminal png \n"
                   <<"set output \"magnTime15.png\"\n"
                   <<"set ylabel \"Magnetism\"\n"
                   <<"set xlabel\"Monte Carlo Step\"\n"
                   <<"plot [][.5:1.5] '-' using 1:2 title \"variation of magnetism for T=1.5\" w l lc 1"<<endl;
                   
magnTime20<<"set title\"Variation of magnetism with respect to time\"\n"
                   <<"set terminal png \n"
                   <<"set output \"magnTime20.png\"\n"
                   <<"set ylabel \"Magnetism\"\n"
                   <<"set xlabel\"Monte Carlo Step\"\n"
                   <<"plot [][.5:1.5] '-' using 1:2 title \"variation of magnetism for T=2.0\" w l lc 1"<<endl;
                   
magnTime225<<"set title\"Variation of magnetism with respect to time\"\n"
                   <<"set terminal png \n"
                   <<"set output \"magnTime225.png\"\n"
                   <<"set ylabel \"Magnetism\"\n"
                   <<"set xlabel\"Monte Carlo Step\"\n"
                   <<"plot [][-.5:1.5] '-' using 1:2 title \"variation of magnetism for T=2.25\" w l lc 1"<<endl;
                   
magnTime40<<"set title\"Variation of magnetism with respect to time\"\n"
                   <<"set terminal png \n"
                   <<"set output \"magnTime40.png\"\n"
                   <<"set ylabel \"Magnetism\"\n"
                   <<"set xlabel\"Monte Carlo Step\"\n"
                   <<"plot [][-.5:.5] '-' using 1:2 title \"variation of magnetism for T=4.0\" w l lc 1"<<endl;
                   
     Corelationdata<<"set title\"Variation of corelation with respect to neighbor\"\n"
                   <<"set terminal png \n"
                   <<"set output \"corelation.png\"\n"
                   <<"set ylabel \"corelation\"\n"
                   <<"set xlabel\"neighbor label\"\n"
                   <<"plot [][] '-' using 1:2 title \"variation of corelation\" w lp ps 1 pt 6"<<endl;
                   
                   
    VarMagnetismPerT<<"set title\"magnetism variance per temperature\"\n"
                   <<"set terminal png \n"
                   <<"set output \"varianceOFmagnetism.png\"\n"
                   <<"set ylabel \"Magnetism variance\"\n"
                   <<"set xlabel\"Temperature\"\n"
                   <<"plot [][-.1:1.2] '-' using 1:2 title \"variance\" w l lc 1"<<endl;

}
/////////////////////////////////////////////////////////////////////////////////////////
// Initialize model, alocate dynamic memory, ...
void Init() {
    cout << "Initialize ..." << endl;
    Randomize();
    s = new int*[L];
    for (int i = 0; i<L; i++) {
      s[i] = new int[L];
      for (int j = 0; j<L; j++)
        s[i][j] = 1;
    }
    
    snapshot << L << endl;
}

// Dealocate dynamic memory, ...
void Done() {
    cout << "Done ..." << endl;

    snapshot.close();

    for (int i = 0; i<L; i++)
      delete[] s[i];
    delete[] s;
}

// Export single snapshot to snapshot.txt
void Export_Single_Snapshot() {
    for (int i = 0; i<L; i++)
      for (int j = 0; j<L; j++)
          if ( s[i][j] > 0 )  
              snapshot << "1" << endl;
          else
              snapshot << "0" << endl;
}

// Rotate i if cross the boundary vs the periodic boundary condition
inline int Period(int i) {
    return i & LMask;
    // or equivalently use: return (i+L) % L;
}

// Calculate Delta E if s(i,j) flipped
float Delta_E_Flip(int i, int j) {
    float sum = 0.0;
    sum = s[i][Period(j+1)] + s[i][Period(j-1)] + s[Period(i-1)][j] + s[Period(i+1)][j];
    return 2*s[i][j]* (h + J * sum);
}

// Perfom single Monte Carlo step
void Single_Monte_Carlo_Step() {
            
        totalMmc=0;
        totalMmc2=0;

    for (int c = 0; c < L*L; c++) { // L^2 infinitesimate step

        int i = Random(L);        // Select random row between [0..L-1]
        int j = Random(L);        // Select random column between [0..L-1]
        float dE = Delta_E_Flip(i, j);
        if ( dE < 0 ) {           // If the total energy decreases after spin-flip, accept the infinitesimal step
            s[i][j] = -s[i][j];
        } else if ( Random() < exp(-dE/T) ) {
            s[i][j] = -s[i][j];
        }
        totalMmc+=s[i][j];
        totalMmc2+=s[i][j]*s[i][j];
    }
    meanMmc=totalMmc/(L*L);
    meanMmc2=totalMmc2/(L*L);
    //cout<<"meanMmc--> "<<meanMmc<<endl;
    // Export_Single_Snapshot();  // Activate this line if you want the snapshot.txt file is being generated during run
}

// Simulate Ising model
void Execute() {
    cout << "Execute the Metropolis algorithm ..." << endl;
            double d=0;
            double Q=0;
   for(T=2;T<3;T+=.01){
       
       cout<<"Temperature --> "<<T<<endl;
       Init();
        for(MCStep=0;MCStep<400;MCStep++){
            Single_Monte_Carlo_Step();
        }
        totalMmc=0;
        totalMmc2=0;
        meanMmc=0;
        meanMmcPerT=0;
        meanMmc2=0;
        meanMmcPerT2=0;

            for(MCStep=0;MCStep<100;MCStep++){
                Single_Monte_Carlo_Step();
                meanMmcPerT+=meanMmc;
                meanMmcPerT2+=meanMmc2;
            }
            meanMmcPerT=meanMmcPerT/100;
            meanMmcPerT2=meanMmcPerT2/100;
            var2M=meanMmcPerT2-meanMmcPerT;
            Sus=var2M/T;
            if(Sus>d){
                d=Sus;
                Q=T;
            }
            varM=sqrt(var2M);
            magnTemp<<T<<"\t"<<meanMmcPerT<<endl;
            VarMagnetismPerT<<T<<"\t"<<varM<<endl;
            Susceptibility<<T<<"\t"<<Sus<<endl;
            
    }
            cout<<"critical T -->"<<Q<<endl;
}

// Perform post statistical analysis
void PostProcess() {
    cout << "Perfom the post statistical analysis ..." << endl;
    // You should perfom the post statistical analysis here
    // .
    // .
    // .
}

// Main routine
int main (int argc, char *argv[]) {
    Init();   
    Initialize_text_file();
    Execute();
    //PostProcess();
    Plot();
    Done();

    cout << "Finish" << endl;
    
    return 0;
}
