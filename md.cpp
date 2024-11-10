#include <cmath>
#include <cstdio>   
#include <cstdlib>  
#include <chrono>   

// 1. constant
const long IMUL = 314159269;
const long IADD = 453806245;
const long MASK = 2147483647;
const double SCALE = 0.4656612873e-9;
static long randSeedP = 17;
constexpr size_t INIT_UCELL_X = 20;
constexpr size_t INIT_UCELL_Y = 20;
constexpr size_t N_MOL = INIT_UCELL_X * INIT_UCELL_Y;  // 400
const int NDIM = 2;          // 2D simulation

// 2. class
class Mol;
class Prop;

// 3. global variable
double region[2];
size_t nMol;
double initUcell[2] = {INIT_UCELL_X, INIT_UCELL_Y};
double velMag;  // velocity magnitude
double vSum[2] = {0.0, 0.0};  // sum of velocities
double virSum;     // virial sum
double uSum;       // potential energy sum
const double epsilon = 1.0;  // LJ parameter epsilon
double rCut;                 // cutoff radius
const double sigma = 1.0;    // LJ parameter sigma
double density; 
double temperature;  
double deltaT;      // time step size
double timeNow;  
size_t stepAvg;   // steps between averaging property measurements
// Add to global variables
size_t stepCount;    // timestep counter

// 4. class
class Mol {
public:
    double r[2];    // location
    double rv[2];   // speed
    double ra[2];   // acceleration
    
    Mol() {
        r[0] = r[1] = 0.0;
        rv[0] = rv[1] = 0.0;
        ra[0] = ra[1] = 0.0;
    }

    Mol(const double r_[2], const double rv_[2], const double ra_[2]) {
        r[0] = r_[0]; r[1] = r_[1];
        rv[0] = rv_[0]; rv[1] = rv_[1];
        ra[0] = ra_[0]; ra[1] = ra_[1];
    }
};

class Prop {
public:
    double val;
    double sum1;
    double sum2;
    
    Prop() : val(0.0), sum1(0.0), sum2(0.0) {}
    Prop(double val_, double sum1_, double sum2_)
        : val(val_), sum1(sum1_), sum2(sum2_) {}
};
void AccumProps(int icode);
void PropZero(Prop &v);
void PropAccum(Prop &v);
void PropAvg(Prop &v, int n);

// 5.global class
Mol mol[N_MOL];
Prop kinEnergy;     // kinetic energy
Prop totEnergy;     // total energy
Prop pressure;      // pressure



inline double Sqr(double x) {
    return x * x;
}

inline double Cube(double x) {
    return x * x * x;
}

//random number generator
double RandR() {
    randSeedP = (randSeedP * IMUL + IADD) & MASK;
    return (randSeedP * SCALE);
}

// random vector
void VRand(double p[2]) {
    double s = 2.0 * M_PI * RandR();
    p[0] = cos(s);
    p[1] = sin(s);
}

// dealing with boundary
void VWrapAll(double v[2]) {
    if (v[0] >= 0.5 * region[0]) {
        v[0] -= region[0];
    } else if (v[0] < -0.5 * region[0]) {
        v[0] += region[0];
    }

    if (v[1] >= 0.5 * region[1]) {
        v[1] -= region[1];
    } else if (v[1] < -0.5 * region[1]) {
        v[1] += region[1];
    }
}

void ApplyBoundaryCond() {
    for (int n = 0; n < nMol; n++) {
        VWrapAll(mol[n].r);
    }
}


void InitCoords() {
    // Calculate gap between molecules
    double gap[2];
    gap[0] = region[0] / initUcell[0];
    gap[1] = region[1] / initUcell[1];
    
    size_t n = 0;
    for (size_t ny = 0; ny < initUcell[1]; ny++) {
        for (size_t nx = 0; nx < initUcell[0]; nx++) {
            // Calculate position: (nx + 0.5, ny + 0.5) * gap - 0.5 * region
            mol[n].r[0] = (nx + 0.5) * gap[0] - 0.5 * region[0];
            mol[n].r[1] = (ny + 0.5) * gap[1] - 0.5 * region[1];
            n++;
        }
    }
}

void InitVels() {
    // Reset velocity sum
    vSum[0] = vSum[1] = 0.0;
    
    // First loop: assign random velocities and compute sum
    for (size_t n = 0; n < nMol; n++) {
        VRand(mol[n].rv);                   
        mol[n].rv[0] *= velMag;           
        mol[n].rv[1] *= velMag;
        vSum[0] += mol[n].rv[0];           
    }
    
    // Second loop: adjust velocities to give zero net momentum
    double vFac = -1.0 / nMol;
    double adjustV[2] = {vSum[0] * vFac, vSum[1] * vFac};
    for (size_t n = 0; n < nMol; n++) {
        mol[n].rv[0] += adjustV[0];
        mol[n].rv[1] += adjustV[1];
    }
}


void InitAccels() {
    for (size_t n = 0; n < nMol; n++) {
        mol[n].ra[0] = mol[n].ra[1] = 0.0;
    }
}


void SetParams() {
    // Calculate cutoff radius
    rCut = pow(2.0, 1.0/6.0 * sigma);
    
    // Define the region size based on density
    region[0] = initUcell[0] / sqrt( density);
    region[1] = initUcell[1] / sqrt( density);
    
    // Set number of molecules
    nMol = N_MOL;
    
    // Calculate velocity magnitude based on temperature
    velMag = sqrt(NDIM * (1.0 - 1.0/nMol) * temperature);

}



void SetupJob() {
    stepCount = 0;
    InitCoords();
    InitVels();
    InitAccels();
    AccumProps(0);   
}

void ComputeForces() {
    double fcVal;           // force value
    double dr[2];          // position difference vector
    double rrCut = Sqr(rCut);
    
    // Reset accelerations and sums
    for (size_t n = 0; n < nMol; n++) {
        mol[n].ra[0] = mol[n].ra[1] = 0.0;
    }
    uSum = virSum = 0.0;
    
    // Compute forces between all pairs
    for (size_t j1 = 0; j1 < nMol - 1; j1++) {
        for (size_t j2 = j1 + 1; j2 < nMol; j2++) {
            // Calculate position difference
            dr[0] = mol[j1].r[0] - mol[j2].r[0];
            dr[1] = mol[j1].r[1] - mol[j2].r[1];
            VWrapAll(dr);
            
            // Calculate squared distance
            double rr = Sqr(dr[0]) + Sqr(dr[1]);
            
            if (rr < rrCut) {
                double r = sqrt(rr);
                
                fcVal = 48.0 * epsilon * pow(sigma, 12) / pow(r, 13) - 
                       24.0 * epsilon * pow(sigma, 6) / pow(r, 7);
                
                // Update accelerations
                mol[j1].ra[0] += fcVal * dr[0];
                mol[j1].ra[1] += fcVal * dr[1];
                mol[j2].ra[0] -= fcVal * dr[0];
                mol[j2].ra[1] -= fcVal * dr[1];
                
                // Update potential energy and virial
                uSum += 4.0 * epsilon * (pow(sigma/r, 12)/r) - pow(sigma/r, 6);
                virSum += fcVal * rr;
            }
        }
    }
}

void LeapfrogStep(int part) {
    if (part == 1) {
        // First part: update velocities by half-step and positions by full step
        for (size_t n = 0; n < nMol; n++) {
            // Update velocity: v(t + Δt/2) = v(t) + (Δt/2)a(t)
            mol[n].rv[0] += 0.5 * deltaT * mol[n].ra[0];
            mol[n].rv[1] += 0.5 * deltaT * mol[n].ra[1];
            
            // Update position: r(t + Δt) = r(t) + Δt*v(t + Δt/2)
            mol[n].r[0] += deltaT * mol[n].rv[0];
            mol[n].r[1] += deltaT * mol[n].rv[1];
        }
    } else {
        // Second part: update velocities by another half-step
        for (size_t n = 0; n < nMol; n++) {
            // Update velocity: v(t + Δt) = v(t + Δt/2) + (Δt/2)a(t + Δt)
            mol[n].rv[0] += 0.5 * deltaT * mol[n].ra[0];
            mol[n].rv[1] += 0.5 * deltaT * mol[n].ra[1];
        }
    }
}

void EvalProps() {
    double vvSum = 0.0;  // sum of velocity squared
    vSum[0] = vSum[1] = 0.0;

    // Calculate velocity sum and velocity squared sum
    for (size_t n = 0; n < nMol; n++) {
        vSum[0] += mol[n].rv[0];
        vSum[1] += mol[n].rv[1];
        vvSum += (mol[n].rv[0] * mol[n].rv[0] + 
                 mol[n].rv[1] * mol[n].rv[1]);
    }

    // Calculate properties
    kinEnergy.val = 0.5 * vvSum / nMol;
    totEnergy.val = kinEnergy.val + (uSum / nMol);
    pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
}

// AccumProps series functions:
void PropZero(Prop &v) {
    v.sum1 = v.sum2 = 0.0;
}

void PropAccum(Prop &v) {
    v.sum1 += v.val;
    v.sum2 += Sqr(v.val);
}

void PropAvg(Prop &v, int n) {
    v.sum1 /= n;
    v.sum2 = sqrt(std::max(v.sum2 / n - Sqr(v.sum1), 0.0));
}

void AccumProps(int icode) {
    switch(icode) {
        case 0:  // Zero properties
            PropZero(totEnergy);
            PropZero(kinEnergy);
            PropZero(pressure);
            break;
            
        case 1:  // Accumulate properties
            PropAccum(totEnergy);
            PropAccum(kinEnergy);
            PropAccum(pressure);
            break;
            
        case 2:  // Average properties
            PropAvg(totEnergy, stepAvg);
            PropAvg(kinEnergy, stepAvg);
            PropAvg(pressure, stepAvg);
            break;
    }
}
void PrintSummary() {
    printf("%zu %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
           stepCount,              // step counter
           timeNow,               // current simulation time
           vSum[0] / nMol,        // average velocity in x direction
           totEnergy.sum1,        // total energy
           totEnergy.sum2,        // total energy standard deviation
           kinEnergy.sum1,        // kinetic energy
           kinEnergy.sum2,        // kinetic energy standard deviation
           pressure.sum1,         // pressure
           pressure.sum2);        // pressure standard deviation
}

void outputMolCoo(const char* filename, size_t mark_1, size_t mark_2) {
    FILE* f = fopen(filename, "w");
    if (!f) {
        printf("Error opening file: %s\n", filename);
        return;
    }

    // Write header
    fprintf(f, "step: %zu\n", stepCount);
    fprintf(f, "ts: %.4f\n", timeNow);
    fprintf(f, "Sigma v: %.4f\n", vSum[0] / nMol);
    fprintf(f, "E: %.4f\n", totEnergy.sum1);
    fprintf(f, "Sigma E: %.4f\n", totEnergy.sum2);
    fprintf(f, "Ek: %.4f\n", kinEnergy.sum1);
    fprintf(f, "Sigma Ek: %.4f\n", kinEnergy.sum2);
    fprintf(f, "P_1: %.4f\n", pressure.sum1);
    fprintf(f, "P_2: %.4f\n", pressure.sum2);
    fprintf(f, "====================\n");

    // Write molecule coordinates
    for (size_t n = 0; n < nMol; n++) {
        if (n == mark_1) {
            fprintf(f, "m-%zu:%.4f,%.4f\n", n, mol[n].r[0], mol[n].r[1]);
        }
        else if (n == mark_2) {
            fprintf(f, "m-%zu:%.4f,%.4f\n", n, mol[n].r[0], mol[n].r[1]);
        }
        else {
            fprintf(f, "o-%zu:%.4f,%.4f\n", n, mol[n].r[0], mol[n].r[1]);
        }
    }
    fprintf(f, "====================");

    fclose(f);
}
void SingleStep() {
    stepCount++;
    timeNow = stepCount * deltaT;

    LeapfrogStep(1);
    ApplyBoundaryCond();
    ComputeForces();
    LeapfrogStep(2);
    EvalProps();
    AccumProps(1);

    if (stepCount % stepAvg == 0) {
        AccumProps(2);
        PrintSummary();
        
        char filename[100];
        sprintf(filename, "coo/%zu.out", stepCount);
        size_t mark_1 = nMol / 2 + nMol / 8;
        size_t mark_2 = mark_1 + 1;
        outputMolCoo(filename, mark_1, mark_2);
        
        AccumProps(0);
    }
}

int main() {
    FILE* infile = fopen("Rap_2_LJP.in", "r");
    if (!infile) {
        printf("Error: cannot open input file\n");
        return 1;
    }
    

    fscanf(infile, "%*s %lf", &deltaT);        // timestep
    fscanf(infile, "%*s %lf", &density);       // density
    fscanf(infile, "%*s %lf", &initUcell[0]);  // initUcell x
    fscanf(infile, "%*s %lf", &initUcell[1]);  // initUcell y
    fscanf(infile, "%*s %zu", &stepAvg);       // step average
    double stepEquil; 
    fscanf(infile, "%*s %lf", &stepEquil);
    double stepLimit;  
    fscanf(infile, "%*s %lf", &stepLimit);
    fscanf(infile, "%*s %lf", &temperature);   // temperature
    fclose(infile);

    printf("steplimit:%.1f\n", stepLimit); 

    system("mkdir -p coo");  
    system("rm -rf coo/*");  
    clock_t start = clock();
    
  
    SetParams();
    SetupJob();
    

    printf("Started MD simulation...\n");
    printf("stepCount timeNow vSum/nMol totE totEsd kinE kinEsd P Psd\n");
    
    bool moreCycles = true;
 
    
    clock_t time_loop_start = clock();
    
    while (moreCycles) {
        SingleStep();
  
        if (stepCount >= (size_t)stepLimit) {
            moreCycles = false;
        }
    }
    
    clock_t time_loop_end = clock();
    
    // 5. 
    double setup_time = (double)(time_loop_start - start) / CLOCKS_PER_SEC;
    double loop_time = (double)(time_loop_end - time_loop_start) / CLOCKS_PER_SEC;
    double total_time = (double)(time_loop_end - start) / CLOCKS_PER_SEC;
    
    printf("\nSimulation completed!\n");
    printf("Setup Time: %.2f seconds\n", setup_time);
    printf("Loop Time: %.2f seconds\n", loop_time);
    printf("Total Time: %.2f seconds\n", total_time);
    
    return 0;
}

