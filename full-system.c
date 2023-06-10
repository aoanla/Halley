#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

//a position vector, units au 
typedef struct {
	double x,y,z;
} vec;

//conversion factor for GM values into au-day-solarmass units
#define GSCALE 2.22972472095e-15

static vec sub(vec a, vec b) {
    vec res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

//modifying additive accumulator
static void accum(vec *a, vec b) {
    a->x += b.x;
    a->y += b.y;
    a->z += b.z;
    return;
}

static double sqr_norm(vec a) {
    return a.x*a.x+a.y*a.y+a.z*a.z;
}

static void mul(vec *a, double b) {
    a->x *= b;
    a->y *= b;
    a->z *= b;
    return;
} 

static vec mul_res(vec a, double b) {
    vec c;
    c.x = a.x*b;
    c.y = a.y*b;
    c.z = a.z*b;
    return c;
}

static void zero(int m, vec a[m]) {
    memset(a, 0, sizeof(vec)*m);
}

//calculate mutual accelerations for a set of objects at positions pons with GM factors Gmasses
//final parameter accns is modified directly to hold the accelerations
static void accns(int num, const vec posns[num], const double Gmasses[num], vec accns[num]) {
    zero(num, accns);
    for (int i=0; i<num; i++) {
        for(int j=i+1; j<num; j++) {
            vec d=sub(posns[i],posns[j]);
            double n = sqr_norm(d); //squared norm
            mul(&d, 1.0/(sqrt(n)*n)); // d / (|d|^2)^3/2 = d_hat / |d|^2
            vec d_ = d;
            mul(&d, 1.0*Gmasses[i]); //by symmetry NOTE I probably have these reversed? d should be i and d_ j?
            accum(accns+j,  d); 
            mul(&d_ , -1.0*Gmasses[j]); //by symmetry
            accum(accns+i, d_);
        }
    } 

    return;    
}                          

#define w0 -1.702414383919315268
#define w1 1.351207191959657634

//Yoshida integrator, 4th-order symplectic 2nd-degree integrator
//integrate positions (posns) and velocities of particles with GM factors Gmasses under acceleration function acc_calc, by a timestep dt
//  posns and velocities are directly updated into the passed arrays
static void integrate(int num, vec posns[num], vec velocities[num], const double Gmasses[num], void (*acc_calc) (int m, const vec[m], const double[m], vec[m]), double dt) {
    static const double c[4] = {w1/2, (w0+w1)/2, (w0+w1)/2, w1/2};
    static const double d[3] = {w1, w0, w1}; 
    vec accns[num];
    for(int iters=0; iters<3; iters++) {
        for(int i=0; i<num;i++) { accum(posns+i, mul_res(velocities[i], c[iters]*dt)); };
        acc_calc(num, posns, Gmasses, accns);
        for(int i=0; i<num;i++) { accum(velocities+i, mul_res(accns[i], d[iters]*dt)); }; 
    }
    for(int i=0; i<num;i++) { accum(posns+i, mul_res(velocities[i], c[3]*dt)); };
    return;
}

//Solve Halley's Comet equation of motion for 200 years after its observation by Edmond Halley in 1682
//print each time it hits perihelion
int main() {

	// to do: these constants
	//			Sol				Jupiter			Saturn			Comet
	vec posns[] = 	{ {2.950832139557744e-3,  -5.425470573959765e-3,  -7.383386124694714e-5},/*Sol*/\
                      {-3.332497263005759,4.112363636591643,5.891307575840340e-2},/*Jupiter*/\
                     {-6.535507530500684,6.380917893848125,1.428260570382428e-1},/*Saturn*/\
			{2.510983241619308e-1,1.950982735843640e-1,-6.755395306466706e-3},/*Mercury*/\
			{-2.620890874732100e-1,-6.808075673769205e-1,6.655856110641488e-3},/*Venus*/\
			{9.630524101520794e-1,-3.128557379506953e-1,-2.422339413846039e-4},/*Earth*/\
			{9.612598229631927e-1,-3.112560037549189e-1,-2.786213887250675e-4},/*Moon*/\
			{-4.415739460116890e-1,1.541057808775840,4.336870492890024e-2},/*Mars*/\
			{1.718669476007500e1,9.933497426989296,-1.866007100527533e-1},/*Uranus*/\
			{2.615297960363933e1,-1.469969306823841e1,-2.995311280685409e-1},/*Neptune*/\
                      {6.391931737411634e-1,-1.121186352054070e-1,1.870482953942138e-1} /*Comet*/};
	vec vels[] = 	{ {	6.838345177814781e-6,5.026348301755031e-6,-2.071993311542051e-7},/*Sol*/\
                      {-5.946882707429968e-3,-4.397581117062436e-3,1.511463721442163e-4},/*Jupiter*/\
                      {-4.208760152990489e-3,-4.000694588560874e-3,2.380878809848479e-4},/*Saturn*/\
			{-2.319387349809163e-2,2.314935923144703e-2,4.029155127823395e-3},/*Mercury*/\
			{1.869504930197275e-2,-7.474624825638375e-3,-1.179108325817677e-3},/*Venus*/\
			{4.975052051056192e-3,1.633653719659161e-2,1.120657612318169e-5},/*Earth*/\
			{4.561144105343879e-3,1.586231559918362e-2,6.830691822388147e-5},/*Moon*/\
			{-1.288189187342671e-2,-2.690962310018411e-3,2.691923720938971e-4},/*Mars*/\
			{-2.000317996107281e-3,3.225767036414139e-3,3.845761974599044e-5},/*Uranus*/\
			{1.517663638954640e-3,2.753873265004750e-3,-9.168067707913533e-5},/*Neptune*/\
                      {-1.551860810762931e-2,-2.496616444169774e-2,4.202270253608422e-4} /*Comet*/	};
	const double Gmasses[] = { GSCALE*132712440041.93938,GSCALE*126686531.900,GSCALE*37931206.234,\
				    GSCALE*22031.86855, GSCALE*324858.592, GSCALE*398600.435436, GSCALE*4902.800066 /**Moon**/,\
				    GSCALE*42828.375214, GSCALE*5793951.256, GSCALE*6835099.97,\
					0.0};

	const int sz = sizeof(Gmasses)/sizeof(Gmasses[0]);
	double  dist, olddist, oldolddist;
	dist = sqr_norm(posns[sz-1]); //sz-1 is the last entry == comet
	olddist = dist;
	oldolddist = olddist;

		
	double dt = 0.1;
	double t = 0;
	const double max_t = 365*200; //200 days
	time_t starttime = mktime(& (struct tm) { .tm_sec = 0, .tm_min = 0, .tm_hour = 0, .tm_mday = 31, .tm_mon = 7, .tm_year = 1682-1900} );
	while(t < max_t) {
		//save old positions
		oldolddist = olddist;
		olddist = dist;

		integrate(sz, posns, vels, Gmasses, accns, dt);
		dist = sqr_norm(posns[sz-1]);
        
		if ( (oldolddist > olddist) && (dist > olddist) ) {
			//perihelion
			//add time conversion
			time_t peritime_t = starttime+(t-dt)*86400;
			struct tm * peritime = localtime(&peritime_t);
			char buffer[100], fmt[55];
			sprintf(fmt, "Perihelion detected at: %%F, time delta %f", t-dt);
			strftime(buffer, 99, fmt, peritime);
			puts(buffer);
		}

		t += dt;
	}
}
