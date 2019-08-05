/*
COMPILE THE CODE WITH: g++ -std=c++11 -g -O3  CH_MC_heatbath.cpp -o MC.o
 (change -mcmodel=large if you want bigger Nh grid)
RUN THE CODE WITH: ./MC.o [spin input file] [temp] [rng_seed]
 ([spin input file]=none if you want to start from a random arrangement)
 
"system.txt" should exist in the directory with information as follows:
j1 j2 jf
u0 umin mmin
eqsweeps avsweeps
 (j1 should be +/-1)
*/
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

//array size for spin lattice
#define n1 20
#define n2 20
#define n3 20
//4th dimension of spin lattice array is 3, to store cartesian vectors
double spins[n1][n2][n3][3] ={};
double J[3]={};
double D = 0;

//grid of points for CDF integral
#define Nh 100001
#define Nm 1024 //ALWAYS 2^N (where N is some integer)
double CDF[Nh][Nm]; //cumulative distribution function

double pi = 3.141592653589793238463;
double integral_bound = 7.0;

//a fudge factor that should be ignored if possible
double x;
//------------------------- Function declarations------------------------------------------
//Monte carlo functions
int PBC(int n, int nmax);
double local_field(double arr[n1][n2][n3][3], double js[3],double D, int i, int j, int k,double h[3]);
double total_X(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin);
double total_U(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin);
double total_energy(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin);


//
double mapping_function(double h[3],double js[3],double hmin, double hmax, double r_theta, double r_phi, double kT,double s_new[3],double s_old[3]);


//Functions for sampling prob dist
double boltz(double h_cos_theta, double m, double u0, double umin, double mmin, double kT); //bolztmann factor
double part_func(double h_cos_theta, double u0, double umin, double mmin, double kT); //1site partition function
void gen_CDF(double hmin, double hmax, double u0, double umin, double mmin, double kT); //does integrals
double int_CDF(double h, double hmin, double hmax, double m); //interpolates CDF
//new functions for the mapping function with an irregular grid
//Non-uniform grid in R space. Points are precisley where they lie on CDF axis. So no information is lost
int R_to_array(int harr, double R);
double int_M_adaptive(double h, double hmin, double hmax, double R);

//Linear algebra helper functions - define my own for now, might import LA library later
double mag(double v[3]);
double dot(double v1[3],double v2[3]);
double cross(double v1[3],double v2[3],double v_out[3]);
double mul(double v1[3],double v2[3],double v_out[3]);
double scalmul(double v[3],double scal);
double matmul(double m1[3][3],double m2[3][3],double m_out[3][3]);
double matvecmul(double m[3][3],double v_in[3]);
int normalise(double v[3]);
double rotate_to_z(double v[3],double m[3][3]);
int transpose(double m[3][3]);
double det(double m[3][3]);
//-----------------------------------------------------------------------------------------------

int main(int argc, char *argv[]){
    //debug
    ofstream debug;
    debug.open("debug.txt");
    
    //command line input arguments
    //cout<<"ah"<<endl;
    string spinin_file=argv[1]; //source of spin structure
    
    double kT=stod(argv[2]);

    double rngseed=stod(argv[3]);
    //seperate mmin used for magnitude distribution only.

    //system parameters to be read from system file
    double Jtot;
    double u0, umin, mmin; //umin is a relic, always zero
    int eqsweeps, avsweeps, sweeps; //equilibriation, averaging and total num of swps
    
    //reading system information
    ifstream systin_str("system.txt");
    systin_str >> J[0] >> J[1] >> J[2] >> D
               >> u0 >> umin >> mmin
               >> eqsweeps >> avsweeps
               >> x;
    
    //Easiest way to modify u0 from bash script
    u0 = stod(argv[4]);
    //D = stod(argv[4]);
    sweeps=eqsweeps+avsweeps;
    //Not sure about this? - only aplicable to triangular lattice?
    //Jtot=2*abs(J[0]) + 2*abs(J[1]) + 6*abs(J[2]);
    //Jtot=2*abs(J[0]) + 2*abs(J[1]) + 2*abs(J[2]);
    Jtot = 6*abs(J[0])+6*abs(D);
    //scaling units for m* in calculation
    //calculation is done in units where m*=1 (historical reasons)
    u0*=(mmin*mmin);
    kT*=(mmin*mmin);
    
    //random number things
    mt19937 rng; //mersenne twister generator
    rng.seed(time(NULL)+100000*rngseed);
    uniform_int_distribution<int> site_picker1(0,n1-1);
    uniform_int_distribution<int> site_picker2(0,n2-1);
    uniform_int_distribution<int> site_picker3(0,n3-1);
    uniform_real_distribution<double> uni_dist(0,1);
    
    //assigning initial spins
    cout << "(*------------*)\n";
    if(spinin_file=="none"){
        cout << "(*) Randomising initial spin config... ";
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                	//randomise cartesian components of spin vectors
                    for(int z=0;z<3;z++){
                    	spins[i][j][k][z] = 2*uni_dist(rng)-1;
                    }
                    //normalise spin vectors
                	//normalise(spins[i][j][k]);	                	
                }
            }
        }
    }
    else{
        ifstream spin_in;
        spin_in.open(spinin_file);
        cout << "(*) Reading initial spin config from "+spinin_file+"... ";
        for(int i=0; i<n1; i++){
            for(int j=0; j<n2; j++){
                for(int k=0; k<n3; k++){
                	for (int z=0;z<3;z++){
                    	spin_in >> spins[i][j][k][z];
                    	//normalise(spins[i][j][k]);
                	}
                }
            }
        }
        spin_in.close();
    }
    cout << "DONE\n";
    

    //setting up probability distribution sampler
    double hmin=0, hmax=integral_bound*Jtot;
    cout << "(*) Generating cumulative distribution function... " << flush;
    gen_CDF(hmin, hmax, u0, umin, mmin, kT);
    cout << "DONE\n";

    //Debug for testing probability distribution
    //ofstream foutp;
    //ofstream foutd ("ha.txt");
    //double i_count = 1000;
    //double s_temp[3] = {1,0,0};
    //double h_vector[3] = {1,0,0};
    //double h_test = 1;
    //for(int i=0; i<=10000; i++){
        //double r = (hmax-hmin)*i/i_count;

    //    double ra = uni_dist(rng);

        //foutd << i << "\t" << h << endl;
        //cout << h << " ";
    //    debug << int_M_adaptive(h_test,hmin,hmax,ra)<<" ";
        //foutp.open("p_sample_"+to_string(i)+".txt");
        //for(int j=0; j<50000; j++){
        //    foutp << (h, hmin, hmax, uni_dist(rng)) << endl;
        //}
        //foutp.close();
    //}
    //cout<<"aaah"<<endl;
    //Initialising averages. In order:
    // energy, squared energy, total spin, total spin squared
    // X energy, squared X energy, U energy, squared U energy
    // total spin on even sites, total spin on odd sites
    // spin (vector) per site, squared (dot product with self) spin per site
    // correlation function (currently out of action)
    double en_avg=0, en2_avg=0, s_avg[3]={}, s2_avg=0;
    double ex_avg=0, ex2_avg=0, eu_avg=0, eu2_avg=0;
    double s_even_avg[3]={}, s_odd_avg[3]={},s2_even_avg=0,s2_odd_avg=0;
    double si_avg[n1][n2][n3][3]={ }, si2_avg[n1][n2][n3]={ };
    double s_corr[n3]={ };
    //Running energy totals
    double enX=total_X(spins, J, u0, umin, mmin);
    double enU=total_U(spins, J, u0, umin, mmin);
    double toten=enX+enU;
    int avsamp=0; //counts the number of data points in average
    
    cout << "(*) Starting Monte Carlo simulation at temperature " << kT/(mmin*mmin) << "... " << flush;

    for(int i=0; i<sweeps; i++){
        for(int j=0; j<n1*n2*n3; j++){
            //pick random site and calculate field
            int a=site_picker1(rng);
            int b=site_picker2(rng);
            int c=site_picker3(rng);
            double s_old[3] ={};
            for(int x=0;x<3;x++){
            	s_old[x]=spins[a][b][c][x]; 
            }
            double mag_old = mag(s_old);
            double h[3] = {};
            local_field(spins, J, D, a, b, c, h);

            //computing initial local energy
            double ex_before=-dot(s_old,h);
            double s_old_sq = dot(s_old,s_old);
            double eu_before= (u0+umin)/(mmin*mmin)*(-2*s_old_sq + s_old_sq*s_old_sq/(mmin*mmin));
            double len_before=ex_before+eu_before;
          
            //pick new spin direction from distribution
            double r_theta=uni_dist(rng);
            double r_phi = uni_dist(rng);
            double s_new[3] = {};
            //Doesn't change magnitude of vectors
            mapping_function(h,J,hmin,hmax,r_theta,r_phi,kT/(mmin*mmin),s_new,s_old);
            //debug<<mapping_function(h,J,hmin,hmax,r_theta,r_phi,kT,s_new)<<" ";
            //debug<<mag(s_new)<<" ";
            
            //pick new magnitude from distribution
            double r = uni_dist(rng);
            double h_scal = abs(dot(h,s_new));
            //debug<<h_scal<<" ";
            double mag_new=int_M_adaptive(h_scal,hmin, hmax, r);
            //cout<<mag_new<<endl;;
            
            scalmul(s_new,mag_new);
            //debug<<mag_new<<endl;



            //compute change in energy and assign new spin
            double ex_after =-dot(s_new,h);
            double s_new_sq = dot(s_new,s_new);
            double eu_after = (u0+umin)/(mmin*mmin)*(-2*s_new_sq + s_new_sq*s_new_sq/(mmin*mmin));
            double len_after= ex_after+eu_after;
            toten+=len_after-len_before;
            enX+=ex_after-ex_before;
            enU+=eu_after-eu_before;
            for(int x=0;x<3;x++){
            	spins[a][b][c][x]=s_new[x];
            }
        }
        //debug << toten << endl;
        if(i>=eqsweeps && i%1==0){
            avsamp++;
            en_avg+=toten;
            ex_avg+=enX;
            eu_avg+=enU;
            en2_avg+=toten*toten;
            ex2_avg+=enX*enX;
            eu2_avg+=enU*enU;
			double s[3]={};
            for(int a=0; a<n1; a++){
                for(int b=0; b<n2; b++){
                    for(int c=0; c<n3; c++){
                    	for(int x=0;x<3;x++){
	                        s[x]=spins[a][b][c][x];
	                        si_avg[a][b][c][x]+=s[x];
                    	}
                        si2_avg[a][b][c]+=dot(s,s);	
//                        if(a==0 && b==0){ //correlation function
//                            for(int m=0; m<n3; m++){
//                                s_corr[c]+=spins[a][b][c]*spins[a][b][PBC(m+c,n3)];
//                            }
//                        }
                    }
                }
            }
        }
        
    }
    cout << "DONE\n(*------------*)\n";
   
    //output streams
    ofstream enavg_out;
    ofstream en2avg_out;
    ofstream eSavg_out;
    ofstream eS2avg_out;
    ofstream savg_out;
    ofstream s2avg_out;
    
    ofstream s_even_avg_out;
    ofstream s_odd_avg_out;
    ofstream s2_even_avg_out;
    ofstream s2_odd_avg_out;

    ofstream siavg_out;
    ofstream si2avg_out;
    ofstream scorr_out;
    ofstream s_out;
    ofstream a_out;
    
    //output files
    enavg_out.open("energy.txt", fstream::app);
    en2avg_out.open("energy2.txt", fstream::app);
    eSavg_out.open("energysplit.txt", fstream::app);
    eS2avg_out.open("energysplit2.txt", fstream::app);
    savg_out.open("spin_total.txt", fstream::app);
    s2avg_out.open("spin2_total.txt", fstream::app);

    s_even_avg_out.open("spin_total_even.txt",fstream::app);
    s_odd_avg_out.open("spin_total_odd.txt",fstream::app);
    s2_even_avg_out.open("spin2_total_even.txt",fstream::app);
    s2_odd_avg_out.open("spin2_total_odd.txt",fstream::app);

    siavg_out.open("spins.txt", fstream::app);
    si2avg_out.open("spins2.txt", fstream::app);
    scorr_out.open("spin_corr.txt", fstream::app);
    s_out.open("spins_after.txt"); //snapshot of spins at the end for resuming
    
    avsweeps=avsamp;
    enavg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << en_avg/avsweeps/(mmin*mmin) << endl;
    en2avg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << en2_avg/avsweeps/(mmin*mmin*mmin*mmin) << endl;
    eSavg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << ex_avg/avsweeps/(mmin*mmin) << "\t" << eu_avg/avsweeps/(mmin*mmin) << endl;
    eS2avg_out << kT/(mmin*mmin) << "\t" << setprecision(12) << ex2_avg/avsweeps/(mmin*mmin) << "\t" << eu2_avg/avsweeps/(mmin*mmin) << endl;    
    siavg_out << kT/(mmin*mmin) << "\t";
    si2avg_out << kT/(mmin*mmin) << "\t";
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                s2_avg+=si2_avg[i][j][k];
                si2avg_out << si2_avg[i][j][k]/avsweeps << "\t";
            	for (int z=0;z<3;z++){
	                s_avg[z]+=si_avg[i][j][k][z];
	                siavg_out << si_avg[i][j][k][z]/avsweeps << "\t";
	                s_out << spins[i][j][k][z] << endl;
	                if(i==0 && j==0){
	                    scorr_out << s_corr[k]/avsweeps/(double)n3 << "\t";
	                }
                }      
                if((i+j+k)%2==0){
                    //even lattice sites
                    s2_even_avg+=si2_avg[i][j][k];
                    for(int z=0;z<3;z++){
                        s_even_avg[z]+=si_avg[i][j][k][z];
                    }
                }else{
                    //odd lattice sites
                    s2_odd_avg+=si2_avg[i][j][k];
                    for(int z=0;z<3;z++){
                        s_odd_avg[z]+=si_avg[i][j][k][z];
                    }
                }
            }
        }
    }
    savg_out << kT/(mmin*mmin) << "\t";
    s_even_avg_out << kT/(mmin*mmin) << "\t";
    s_odd_avg_out << kT/(mmin*mmin) << "\t";

    for(int x=0;x<3;x++){
    	s_avg[x]/=(avsweeps);
        s_even_avg[x]/=(avsweeps);
        s_odd_avg[x]/=(avsweeps);
    	savg_out << " " << s_avg[x];
        s_even_avg_out << " " << s_even_avg[x];
        s_odd_avg_out << " " << s_odd_avg[x];
    }
    savg_out << endl;
    s_even_avg_out<<endl;
    s_odd_avg_out<<endl;

    //s2avg_out << kT/(mmin*mmin) << "\t" << s2_avg/avsweeps << endl;
    //output both average of spin magnitudes^2 and magnitude^2 of average spin
    s2avg_out << kT/(mmin*mmin) << "\t" << s2_avg/(avsweeps) <<" "<< dot(s_avg,s_avg) << endl;
    s2_even_avg_out << kT/(mmin*mmin) << "\t" << s2_even_avg/(avsweeps) <<" "<< dot(s_even_avg,s_even_avg) << endl;
    s2_odd_avg_out << kT/(mmin*mmin) << "\t" << s2_odd_avg/(avsweeps) <<" "<< dot(s_odd_avg,s_odd_avg) << endl;

    enavg_out.close();
    en2avg_out.close();
    savg_out.close();
    s2avg_out.close();
    s_even_avg_out.close();
    s_odd_avg_out.close();
    s2_even_avg_out.close();
    s2_odd_avg_out.close();
    siavg_out.close();
    si2avg_out.close();
    scorr_out.close();
    s_out.close();
    debug.close();

    //return 0;
}

//***********************************//
//-----------------------------------//
//       Monte Carlo FUNCTIONS       //
//-----------------------------------//
//***********************************//

//NEED TO EDIT STUFF BELOW HERE

//PBC=Periodic Boundary Conditions
int PBC(int n, int nmax){
    while(n>=nmax){ n-=nmax; }
    while(n<0){ n+=nmax; };
    return n;
}

//Computes local field. This is where lattice information enters
double local_field(double arr[n1][n2][n3][3], double js[3],double D, int i, int j, int k,double h[3]){
    i=PBC(i,n1);
    j=PBC(j,n2);
    k=PBC(k,n3);

    //Cubic lattice
    for (int x=0;x<3;x++){
    	h[x]=js[0]*(arr[i][j][PBC(k+1,n3)][x]+arr[i][j][PBC(k-1,n3)][x]) +
    		 js[1]*(arr[i][PBC(j+1,n2)][k][x]+arr[i][PBC(j-1,n2)][k][x]) +
    		 js[2]*(arr[PBC(i+1,n1)][j][k][x]+arr[PBC(i-1,n1)][j][k][x]);
    }
    

    //Edit here to add DM interactions
    double dms[6][3] = {};
    //n_ij is the set of unit vectors from a site to any neighbouring site. Trivial for cubic but
    //likely more complicated for triangular lattice
    double n_ij[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    cross(arr[i][j][PBC(k+1,n3)],n_ij[0],dms[0]);
    cross(arr[i][j][PBC(k-1,n3)],n_ij[1],dms[1]);
    cross(arr[i][PBC(j+1,n2)][k],n_ij[2],dms[2]);
    cross(arr[i][PBC(j-1,n2)][k],n_ij[3],dms[3]);
    cross(arr[PBC(i+1,n1)][j][k],n_ij[4],dms[4]);
    cross(arr[PBC(i-1,n1)][j][k],n_ij[5],dms[5]);

    for (int y=0;y<6;y++){
        for (int x=0;x<3;x++){
            h[x]+=D*dms[y][x];
        }
    }


    //include next nearest neighbours
    for (int x=0;x<3;x++){
        h[x]-=0.0625*(js[0]*(arr[i][j][PBC(k+2,n3)][x]+arr[i][j][PBC(k-2,n3)][x]) +
                      js[1]*(arr[i][PBC(j+2,n2)][k][x]+arr[i][PBC(j-2,n2)][k][x]) +
                      js[2]*(arr[PBC(i+2,n1)][j][k][x]+arr[PBC(i-2,n1)][j][k][x]));
    }
    cross(arr[i][j][PBC(k+2,n3)],n_ij[0],dms[0]);
    cross(arr[i][j][PBC(k-2,n3)],n_ij[1],dms[1]);
    cross(arr[i][PBC(j+2,n2)][k],n_ij[2],dms[2]);
    cross(arr[i][PBC(j-2,n2)][k],n_ij[3],dms[3]);
    cross(arr[PBC(i+2,n1)][j][k],n_ij[4],dms[4]);
    cross(arr[PBC(i-2,n1)][j][k],n_ij[5],dms[5]);

    for (int y=0;y<6;y++){
        for (int x=0;x<3;x++){
            h[x]-=0.125*D*dms[y][x];
        }
    }

    //Triangle lattice
    /*
    for (int x=0;x<3;x++){
    	h[x] = (js[0]*(arr[i][j][PBC(k+1,n3)][x] + arr[i][j][PBC(k-1,n3)][x]) +
            	js[1]*(arr[i][j][PBC(k+2,n3)][x] + arr[i][j][PBC(k-2,n3)][x]) +
            	js[2]*(arr[PBC(i+1,n1)][j][k][x] + arr[PBC(i-1,n1)][j][k][x] +
                   arr[i][PBC(j+1,n2)][k][x] + arr[i][PBC(j-1,n2)][k][x] +
                   arr[PBC(i+1,n1)][PBC(j-1,n2)][k][x] + arr[PBC(i-1,n1)][PBC(j+1,n2)][k][x]));
    }
    */
    //return 0;
}

//Computes TOTAL energy, computed using local_field function
double total_energy(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin){
    double e=0;
    double m[3] ={};
    double l_field[3] = {};
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
                for(int x=0;x<3;x++){
                	m[x]=arr[i][j][k][x];
                }
                local_field(arr, js, D, i, j, k,l_field);
                double m_sq = dot(m,m);
                e+=u0+(-0.5*dot(m,l_field) + (u0+umin)/(mmin*mmin)*(-2*m_sq + m_sq*m_sq/(mmin*mmin)));
            }
        }
    }
    return e;
}

//Computes EXCHANGE energy
double total_X(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin){
    double ex=0;
    double m[3]={};
	double l_field[3] = {};

    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
            	for(int x=0;x<3;x++){
                	m[x]=arr[i][j][k][x];
            	}
                local_field(arr, js,D, i, j, k,l_field);
                ex+=-0.5*dot(m,l_field);
            }
        }
    }
    
    return ex;
}

//Compute ON-SITE energy
double total_U(double arr[n1][n2][n3][3], double js[3], double u0, double umin, double mmin){
    double eu=0;
    double m[3]={0,0,0};
	double l_field[3] = {0,0,0};
    for(int i=0; i<n1; i++){
        for(int j=0; j<n2; j++){
            for(int k=0; k<n3; k++){
            	for(int x=0;x<3;x++){
                	m[x]=arr[i][j][k][x];
                }
                local_field(arr, js,D, i, j, k,l_field);
                double m_sq = dot(m,m);
                eu+=u0+((u0+umin)/(mmin*mmin)*(-2*m_sq + m_sq*m_sq/(mmin*mmin)));
            }
        }
    }
    
    return eu;
}


double mapping_function(double h[3], double js[3],double hmin, double hmax, double r_theta, double r_phi, double kT,double s_new[3],double s_old[3]){
	//Analytic mapping function for inverse transform sampling
	//Rotate h to be along z axis, and store rotation matrix. Angles are easy to define wrt z axis.
	//Generate random vector from angular distribution and rotate back.
	
	//initialise s_new output at 0
	for(int x=0;x<3;x++){
		s_new[x]=0;
	}
    double old_mag = mag(s_old);
	double h_mag = mag(h)*old_mag;
	//rotation maps e_z to direction of h. 
	double rotation[3][3] = {};
	rotate_to_z(h,rotation);
    //cos(theta) component. Don't bother with arccos as it's only going to be used with trig funcs
	double cos_theta = (kT/(h_mag))*log((1-r_theta)*exp(h_mag/kT)+r_theta*exp(-h_mag/kT));
	double sin_theta = sqrt(1-cos_theta*cos_theta); 

    //error catching for very low temperatures
    //as T -> 0, sometimes cos_theta diverges
	if(std::isnan(cos_theta)){
		cos_theta=1;
		sin_theta=0;
	}

	//phi component
	double phi = 2*((double)pi*r_phi);
	//Set s_new to vector with random angles drawn from specified distribution. 
	//Angles wrt e_z (i.e. typical spherical polar coords)

	s_new[0] = cos(phi)*sin_theta;
	s_new[1] = sin(phi)*sin_theta;
	s_new[2] = cos_theta;
	matvecmul(rotation,s_new);
    
	//return det(rotation);
}


//***********************************//
//-----------------------------------//
// Linear Algebra helper functions   //
//-----------------------------------//
//***********************************//


double mag(double v[3]){
	return (sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

double dot(double v1[3],double v2[3]){
	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

double cross(double v1[3],double v2[3],double v_out[3]){
    v_out[0]=v1[1]*v2[2]-v1[2]*v2[1];
    v_out[1]=v1[2]*v2[0]-v1[0]*v2[2];
    v_out[2]=v1[0]*v2[1]-v1[1]*v2[0];  
}

double mul(double v1[3],double v2[3],double v_out[3]){
    for(int x=0;x<3;x++){
        v_out[x]=v1[x]*v2[x];
    }
}

double scalmul(double v[3],double scal){
    //in place scalar multiplication of vector
    for(int x=0;x<3;x++){
        v[x]*=scal;
    }
}

double matvecmul(double m[3][3],double v_in[3]){
  //multiply vector in place by matrix
  double v_out[3] = {0,0,0};
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++){
      v_out[x]+=m[x][y]*v_in[y];
    }
  }
  for(int x=0;x<3;x++){
    v_in[x]=v_out[x];
  }
}

double matmul(double m1[3][3],double m2[3][3],double m_out[3][3]){
  //naive matrix multiplication
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++){
      m_out[x][y] = 0;
    }
  }
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++){
      for(int z=0;z<3;z++){
        m_out[x][y]+=m1[x][z]*m2[z][y];
      }
    }
  }
}

int normalise(double v[3]){
  //normalise vector in place
  double mag = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  for(int x=0;x<3;x++){
    v[x] = v[x]/mag;
  }

}

int transpose(double m[3][3]){
  //transpose matrix in place
  double m_copy[3][3] = {};
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++){
      m_copy[x][y] = m[x][y];
    }
  }
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++){
      m[x][y] = m_copy[y][x];
    }
  }
}

double det(double m[3][3]){
	return (m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
		   -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
		   +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]));
}



double rotate_to_z(double v[3],double m[3][3]){
  //Calculates rotation matrix required to tranform e_z to given vector v.
  //First find rotation matrix mapping v to e_z, then invert(transpose).

  //initialise stuff
  double vec[3] = {};
  for(int x=0;x<3;x++){
    vec[x]=v[x];
    for(int y=0;y<3;y++){
      m[x][y] = 0;
    }
  }
  normalise(vec);
  double x = vec[0];
  double y = vec[1];
  double z = vec[2];

  //If y is zero, there is divide by 0 error
  if (y!=0){
    double k = x/y;
    double k_root = sqrt(k*k+1);
    double rz[3][3] = {{1/k_root,-k/k_root,0},{k/k_root,1/k_root,0},{0,0,1}};
    double rx[3][3] = {{1,0,0},{0,z,-copysign(sqrt(1-z*z),y)},{0,copysign(sqrt(1-z*z),y),z}};
    //Chain together rotation about z and x axes.
    matmul(rx,rz,m);
  }
  else{
  	//Only need rotation about y axis
    double ry[3][3] = {{z,0,-copysign(sqrt(1-z*z),x)},{0,1,0},{copysign(sqrt(1-z*z),x),0,z}};
    for(int a=0;a<3;a++){
      for(int b=0;b<3;b++){
        m[a][b] = ry[a][b];
      }
    }
  }
  //Transpose matrix as inverse rotation (i.e. e_z -> v) is needed.
  transpose(m);
}








//***********************************//
//-----------------------------------//
//        M sampler FUNCTIONS        //
//-----------------------------------//
//***********************************//

//boltzmann factor
double boltz(double h_cos_theta, double m, double u0, double umin, double mmin, double kT){
    //double m_mag = mag(m);
    //return exp( -(-h_cos_theta*m +x+ u0+(u0+umin)*m*m/(mmin*mmin)*(m*m/(mmin*mmin)-2))/kT );
    //return m*m*exp( -(x+u0-h_cos_theta*m+(u0+umin)*m*m/(mmin*mmin)*(m*m/(mmin*mmin)-2))/kT);
    return m*m*exp( -(x+u0+m*(-h_cos_theta + (u0+umin)*(-2*m + m*m*m/(mmin*mmin))/(mmin*mmin)))/kT );
}

//partition function integral over m
// uses trapezium method
double part_func(double h_cos_theta, double u0, double umin, double mmin, double kT){
    //change this?
    //double dm=2.0/(Nm-1);
    double dm=integral_bound/(Nm-1);
    double Zmid=0;
    //Need to modify integrand to account for spherical integration
    for(int i=0; i<Nm; i++){
        Zmid+=boltz(h_cos_theta, i*dm, u0, umin, mmin, kT);
    }
    double z=0.5*dm*(boltz(h_cos_theta, 0, u0, umin, mmin, kT)+2*Zmid+boltz(h_cos_theta, integral_bound, u0, umin, mmin, kT));
    return z;
}

//integrates probability (boltz/Z) to get cumulative dist function
void gen_CDF(double hmin, double hmax, double u0, double umin, double mmin, double kT){
    double h, dh=(hmax-hmin)/(Nh-1); //because we need to data points at either end!
    //double dm=2.0/(Nm-1);
    double dm =integral_bound/(Nm-1);
    double Z, cdf;
    
    for (int i=0; i<Nh; i++){
        h=hmin+i*dh;
        //calculate partition function for this h
        Z=part_func(h, u0, umin, mmin, kT);
        //calculate CDF
        // Lots of safety features in here, best to leave them be
        cdf=0;
        CDF[i][0]=0; //force the first point to be 0
        for(int j=1; j<Nm; j++){
            //Integrating from 0 to 1
            cdf+=0.5*dm*(boltz(h, (j-1)*dm, u0, umin, mmin, kT)
                         + boltz(h, j*dm, u0, umin, mmin, kT));
            if(cdf==0){
                CDF[i][j]=0;
            }
            else if(std::isinf(cdf)==1){
                CDF[i][j]=1;
            }
            else{
                if(std::isnan(Z)==1){
                    cout << "NOT DONE\n Z is nan, pick better X\n";
                    exit(-1);
                }
                else if(std::isinf(Z)==1){
                    cout << "NOT DONE\n Z is infinite, pick better X\n";
                    exit(-1);
                }
                else if(Z==0){
                    cout << "NOT DONE\n Z is zero, pick better X\n";
                    exit(-1);
                }
                CDF[i][j]=cdf/Z;
            }
        }
        for(int j=0; j<Nm; j++){
            CDF[i][j]/=CDF[i][Nm-1]; //force the last point to be 1 by normalising
            //cout<<CDF[i][j]<<" ";
        }
        //cout<<endl;
        //cout<<CDF[i][20]<<" ";
    }
    
}

//Bilinear interpolation of the CDF array
//RELIC
double int_CDF(double h, double hmin, double hmax, double m){
    // find grid points either side of h coord
    double dh=(hmax-hmin)/(Nh-1);
    double harr=(h-hmin)/dh;
    int h1arr=floor(harr), h2arr=ceil(harr);
    
    //find grid points either side of m coord
    double dm=integral_bound/(Nm-1);
    double marr=(m+1.0)/dm;
    int m1arr=floor(marr), m2arr=ceil(marr);
    
    double outp;
    
    if(h1arr!=h2arr){
        if(m1arr!=m2arr){
            double m1=-1+dm*m1arr, m2=-1+dm*m2arr;
            double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
            outp = (CDF[h1arr][m1arr]*(h2-h)*(m2-m) +
                    CDF[h2arr][m1arr]*(h-h1)*(m2-m) +
                    CDF[h1arr][m2arr]*(h2-h)*(m-m1) +
                    CDF[h2arr][m2arr]*(h-h1)*(m-m1))/(h2-h1)/(m2-m1);
        }
        else if(m1arr==m2arr){
            double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
            outp = CDF[h1arr][(int)marr] + (CDF[h2arr][(int)marr]-CDF[h1arr][(int)marr])*(h-h1)/(h2-h1);
        }
    }
    else{
        if(m1arr!=m2arr){
            double m1=-1+dm*m1arr, m2=-1+dm*m2arr;
            outp = CDF[(int)harr][m1arr] + (CDF[(int)harr][m2arr]-CDF[(int)harr][m1arr])*(m-m1)/(m2-m1);
        }
        else if(m1arr==m2arr){
            outp = CDF[(int)harr][(int)marr];
        }
    }
    
    return outp;
}

//Converts a given value of R into array units. Gives the array value below it.
// Given a random number 0<R<1, this scans through the CDF and finds the
// grid point in the CDF array below where this R corresponds too
int R_to_array(int harr, double R){
    int lowi=0, highi=Nm-1;
    int midhighi=(highi+1)/2, midlowi=midhighi-1;
    double lowR, midlowR, midhighR, highR;
    //Uses a bisection type algorithm, which is why
    // Nm should always be a power of 2
    while(highi-lowi!=1){
        lowR=CDF[harr][lowi];
        highR=CDF[harr][highi];
        midlowR=CDF[harr][midlowi];
        midhighR=CDF[harr][midhighi];
        
        if(R>lowR && R<midlowR){ //if it's in the bottom half
            highi=midlowi;
            midhighi=lowi+(highi+1-lowi)/2;
            midlowi=midhighi-1;
        }
        else if(R>midhighR && R<highR){ //if it's in the top half
            lowi=midhighi;
            midhighi=lowi+(highi+1-lowi)/2;
            midlowi=midhighi-1;
        }
        else if(R>midlowR && R<midhighR){ //if it's in the middle
            lowi=midlowi;
            highi=midhighi;
        }
    }
    
    return lowi;
}

//Mapping function. Returns m for any R and h by linearly interpolating array
// This inverts the CDF *on the fly* using R_to_array. Four points to be
// interpolated between are now *non-rectangular*, so need to use something fancy
double int_M_adaptive(double h, double hmin, double hmax, double R){
    //magnitutude of s_new should be 1, as only a direction has been assigned so far
    //cout<<h<<" ";
    double dm=integral_bound/(Nm-1.0);
    double dh=(hmax-hmin)/(Nh-1);
    double harr=(h-hmin)/dh;
    int h1arr=floor(harr), h2arr=ceil(harr);
    
    double outp;
    
    if(h1arr!=h2arr){
        //converting R into array units
        int R11arr=R_to_array(h1arr, R);
        int R12arr=R11arr+1;
        int R21arr=R_to_array(h2arr, R);
        int R22arr=R21arr+1;
        //all the coords we need
        double h1=hmin+dh*h1arr, h2=hmin+dh*h2arr;
        double R11=CDF[h1arr][R11arr], R12=CDF[h1arr][R12arr];
        double R21=CDF[h2arr][R21arr], R22=CDF[h2arr][R22arr];
        
        //now do some algebra with ^^those^^ things to interpolate
        //convenient definitions
        double Dh=h2-h1;
        double DR1=R12-R11, DR2=R22-R21;
        //quadratic coefficients
        double a=Dh*(DR1-DR2);
        double b=(h1-h)*(DR1-DR2)-Dh*DR1;
        double g=(h-h1)*DR1;
        
        //finding s and t
        double s, s1, s2, t;
        if(a==0){ s=-g/b;} //don't need to solve quadratic (often true to double prec)
        else{ //need to solve quadratic
            double ss=sqrt(b*b-4*a*g);
            s1=(-b+ss)/(2*a), s2=(-b-ss)/(2*a);
            if(s1>0 && s1<1){ s=s1; }
            else{ s=s2; }
        }
        t=(R-R12-s*(R22-R12))/(R11-R12+s*(DR1-DR2));
        
        //bilinear interpolation
        //outp = (-1+dm*R12arr)*(1-s)*(1-t) +
        //(-1+dm*R22arr)*s*(1-t) +
        //(-1+dm*R11arr)*(1-s)*t +
        //(-1+dm*R21arr)*s*t;
        outp = (dm*R12arr)*(1-s)*(1-t) +
        (dm*R22arr)*s*(1-t) +
        (dm*R11arr)*(1-s)*t +
        (dm*R21arr)*s*t;
    }
    else{//(if we lie directly on a grid point) NEVER EVEN GETS USED SO IS INCOMPLETE
        int R1arr=R_to_array(h1arr, R);
        int R2arr=R1arr+1;
        
        double M1=dm*R1arr;
        double M2=dm*R2arr;
        double R1=CDF[h1arr][R1arr];
        double R2=CDF[h1arr][R2arr];
        outp = M1 + (R-R1)*(M2-M1)/(R2-R1);
    }
    
    return outp;
}


