
/******************************************************************************/
/*** Computational Physics Course                                           ***/
/*** Name :   Mohammad Hemati                                               ***/
/*** Stu Num:        951133                                                 ***/
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
#include <sstream>

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
double T=2.48;                         // Tempreture
ofstream snapshot("snapshot.txt", // Keeps model snapshot
                  ios::out | ios::trunc);  
////////////////////////////////////////////////////////////////////////////////////////
//global file textt

ofstream Corelationdata248("Corelationdata248.txt");

////////////////////////////////////////////////////////////////////////////////////////
//define some global variables
double totalMmc=0;              //total magnetism in each Monte-Carlo step
double meanMmc=0;               //mean of magnetism in each Monte-Carlo step
int MCStep=0;               
double meanMmcPerT=0;           //mean of magnetism
double meanMmc2=0;              //mean of magnetism^2 in each MC step
double totalMmc2=0;             //total magnetism*magnetism in each MC step
double meanMmcPerT2=0;          //mean of magnetism*magnetism 
double varM=0;                  //variance of magnetism
double var2M=0;                 //variance of m^2
double Sus=0;                   //susceptibility of system
const int Sweep=1000;           //total MC steps that sweeps for get physical data
double Corl15[128];             //corelation for spin128 in Tempreture 1.5
double Corl225[128];            //corelation for spin128 in Tempreture 2.25
double Corl248[128];             //corelation for spin128 in Tempreture 2.5
double Corl30[128];             //corelation for spin128 in Tempreture 3.0
double Corl20[128];             //corelation for spin128 in Tempreture 2.0

////////////////////////////////////////////////////////////////////////////////////////
//prototype function
void Initialize_text_file();      //create text file belong whit some gnuplot commands
void Corelation();                //Calculate the corelation of system 
void Randomize_Spins();           //initial direction of system gonna random
void Plot();                      //Call gnuplot by system command
void Init();                      // Initialize model, alocate dynamic memory, ...
void Execute();                   // Simulate Ising model
void PostProcess();               // Perform post statistical analysis
void Done();                      // Dealocate dynamic memory, ...
/////////////////////////////////////////////////////////////////////////////////////////
//........................My Interface.................................................//
//.....................................................................................//
//Plot function
void Plot(){
    cout<<"Plotting ....."<<endl;
    Corelationdata248.close();
        cout<<"Corelationdata248"<<endl;
    system("gnuplot 'Corelationdata248.txt'");

}

//OUTput file
void Initialize_text_file(){

                   
                   
    Corelationdata248<<"set title\"Variation of Corelation with respect to neighbors\"\n"
                   <<"set terminal png \n"
                   <<"set output \"Corelationdata248.png\"\n"
                   <<"set ylabel \"corelation\"\n"
                   <<"set xlabel\"number of neighbor\"\n"
                   <<"plot [][] '-' using 1:2 title \"variation of corelation for T=2.48 Critical T\" w l lc 1"<<endl;                   
                   


}

//Give random value to spins
void Randomize_Spins(){
    for (int i = 0 ; i < L ; i++) 
    {
        for (int j = 0 ; j < L ; j++)  
        {
            s[i][j] = (rand() % 2) * 2 - 1;
        }
    }
}

//Corelation function
void Corelation(){
    int c=128;
    for(int i=1;i<(L/2);i++){
            Corl248[i]+=s[c][c]*(s[c-i][c]+ s[c+i][c]+ s[c][c-i]+ s[c][c+i]);
    }
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
    for (int i = 0; i<L; i++, snapshot << endl)
      for (int j = 0; j<L; j++)
          if ( s[i][j] > 0 )  
              snapshot << "1" << "\t";
          else
              snapshot << "0" << "\t";
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
        
    }
    for(int k=0; k<L;k++){
        for(int p=0; p<L;p++){
            totalMmc+=s[k][p];
        }
    }
    meanMmc=totalMmc/(L*L);

}

// Simulate Ising model 
void Execute() {
    cout << "Execute the Metropolis algorithm ..." << endl;
    
       cout<<"Temperature --> "<<T<<endl;
       Init();
        for(MCStep=0;MCStep<400;MCStep++){                          //this loop is for equlibirum
            Single_Monte_Carlo_Step();
        }
        
            for(MCStep=0; MCStep<Sweep; MCStep++){                  //this loop is run for giving physical values
                Single_Monte_Carlo_Step();

                    Corelation();                                   //call corelation function 
            }
            
            //store corelation data for seperate temperature
            for(int i=0; i<L/2; i++){
                    Corelationdata248<<i<<"\t"<<Corl248[i]/(4*Sweep)<<endl;
            }
}

// Perform post statistical analysis
void PostProcess() {
    cout << "Perfom the post statistical analysis ..." << endl;
    
}

// Main routine
int main (int argc, char *argv[]) {
    Init(); 
    //Randomize_Spins();
    Initialize_text_file();
    Execute();
    //PostProcess();
    Plot();
    Done();
    cout << "Finish" << endl;
    
    return 0;
}
