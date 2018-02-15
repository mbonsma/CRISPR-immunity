#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <random>
#include <cmath> 
#include <Eigen/Dense>

//Code to simulate the bacteria-phage system with CRISPR.

// Generates a uniform distribution of numbers between [0,1[ (double). To pull a number from distribution, use ud(mrand)
std::mt19937 mrand(std::time(0));  // seed
std::uniform_real_distribution<> ud(0,1);

unsigned long  poisson(double lambda)
{ //Function generates a poisson distribution with mean lambda and returns one random integer from that distribution.
    std::poisson_distribution<unsigned long > d(lambda);
    return d(mrand);
}

std::string filename(const int m, const float xi, const float eta, const float R,  const float pv0, const int exp, const int c0Exp, const int repeat, const float xinit, const float yinit, const float zinit, const float nuinit)
{  //Function to generate the name for the output text file. 
	//std::string y = "./data/";
	std::string y="./data/";
	y += std::to_string(m);
	y += ("-");
	y += std::to_string(xi);
	y += ("-");
	y += std::to_string(eta);
	y += ("-");
	y += std::to_string(R);
	y += ("-");
	y += std::to_string(pv0);
	y += ("-");
	y += std::to_string(exp);
	y += ("-");
	y += std::to_string(c0Exp);
	y += ("-");
	y += std::to_string(repeat);
	y += ("-");
	y += std::to_string(xinit);
	y += ("-");
	y += std::to_string(yinit);
	y += ("-");
	y += std::to_string(zinit);
	y += ("-");
	y += std::to_string(nuinit);
	y += (".txt");
	return y;
}

long int jumpF(const int exp, const int maxDataPoints)
{//Function to calculate how many iterations to skip when saving data. There is a maximum of 10^maxDataPoints iterations saved, and a maximum of 10^exp iterations. 
	if((exp-maxDataPoints)>=1){
		return pow(10,exp-maxDataPoints);
	}
	else {
		return 1; //default minimum
  	}
}

int NF(const int m)
{//Number of species (nbo,nbi,Nv,C)	
	return m+3;
}

int MF(const int m)
{//Number of reactions
	return 6*m+7;
}

void changeToPopulation(const int m, const int B, Eigen::MatrixXi& v)
{   //change in population. v(i,j) gives the change in species i due to reaction j.
	v.setZero();

    for(int i=0;i<m+1;++i){
        //bacteria growth
        v(i,i)=1; 
        v(m+2,i)=-1;
        //bacteria flow
        v(i,i+m+1)=-1;
 
    }
    //phage flow
    v(m+1,2*m+2)=-1;
    //carbon in
    v(m+2,2*m+3)=1;
    //carbon out
    v(m+2,2*m+4)=-1;
 
    for(int i=0;i<m;++i){
        //loose spacer
        v(i+1,2*m+5+i)=-1;
        v(0,2*m+5+i)=1;
        //v-b0 gain spacer
        v(0,3*m+6+i)=-1;
        v(m+1,3*m+6+i)=-1;
        v(i+1,3*m+6+i)=1;
        //v-bi, b wins
        v(m+1,4*m+7+i)=-1;
        //v-bi, v wins
        v(i+1,5*m+7+i)=-1;
        v(m+1,5*m+7+i)=B-1;
    }
    //interaction v-b0, b wins
    v(m+1,3*m+5)=-1;
    //interaction v-b0, v wins
    v(0,4*m+6)=-1;
    v(m+1,4*m+6)=B-1;
}

void initialPopulation(Eigen::VectorXd& x, const double nbi, const double nbsi, const double nvi, const double ci, const int m)
{//Initial population vector.
	x.setZero();
	x(0)=nbi;
	x(1)=nbsi;
	x(m+1)=nvi;
	x(m+2)=ci;	
}

void totalRates(const Eigen::VectorXd& x, Eigen::VectorXd& a, const int m, const float g, const float f, const double c0, const int B, const float r, const float alpha, const float pv0, const float eta, const float e)
{//Calculates the total rates of the reactions a(j) for the current iteration.
	for(int i=0;i<m+1;++i){
		//bacteria growth
		a(i)=g*x(m+2)*x(i); 
		//bacteria flow
		a(i+m+1)=f*x(i);
	}
	//phage flow
	a(2*m+2)=f*x(m+1);
	//carbon in
	a(2*m+3)=f*c0;
	//carbon out
	a(2*m+4)=f*x(m+2);

	for(int i=0;i<m;++i){
		//loose spacer
		a(2*m+5+i)=r*x(i+1);
		//v-b0 gain spacer
		a(3*m+6+i)=alpha*(1-pv0)*eta/m*x(0)*x(m+1);
		//v-bi, b wins
		a(4*m+7+i)=alpha*(1-pv0+pv0*e)*x(i+1)*x(m+1);
		//v-bi, v wins
		a(5*m+7+i)=alpha*pv0*(1-e)*x(i+1)*x(m+1);
	}
	//interaction v-b0, b wins, no spacer
	a(3*m+5)=alpha*(1-pv0)*(1-eta)*x(0)*x(m+1);
	//interaction v-b0, v wins
	a(4*m+6)=alpha*pv0*x(0)*x(m+1);
}

void saveToFile(const Eigen::VectorXd& x, std::ofstream& file, const int m, const int n, const double t, const int method)
{//Function to save the data from the current iteration to the file.
	int N=NF(m);
	file << n << "\t";
	for(int i=0;i<N;++i){
		file << std::scientific << x(i) << "\t" ;
	}
	file << std::scientific << t << "\t";
	file <<   method << "\n";
} 


    
double totalBacteria(const Eigen::VectorXd& x, const int m)
{//Total number of bacteria with and without spacer.
	double Nb=0;
	for(int i=0;i<m+1;++i){
		Nb += x(i);
	}
	return Nb;
}

double tauTimeStep(const Eigen::MatrixXi& v, const Eigen::VectorXd& a, const Eigen::VectorXd& x, const int m, const float epsilon)
{//Calculates tau, the time interval to jump over. It's the largest time step such that the rates of reactions don't change too much.
	const int N=NF(m);
	const int M=MF(m);
	std::vector<double> taup; //taup stores all possible values of tau
	Eigen::VectorXd mu(N);
	Eigen::VectorXd sig2(N);
	mu.setZero();
	sig2.setZero();
	taup.reserve(2*N);
	double d1;
	const double d2=1.;

	for (int i=0; i<N ;++i){ 
		for(int j=0;j<M;++j){
			mu(i) += v(i,j)*a(j);
			sig2(i) += v(i,j)*v(i,j)*a(j);
		}
		d1 =epsilon*x(i)/2.; 
		taup.push_back(std::max(d1,d2)/std::abs(mu(i))); //push_back adds an element at the end of a vector
		taup.push_back(std::pow(std::max(d1,d2),2)/sig2(i));
		
	}
	std::vector<double>::iterator min = std::min_element(std::begin(taup), std::end(taup)); 
	return *min;
}

double gillespieTime(double a0)
{//calculates time until next reaction.
	double random=ud(mrand);
	return 1./a0*log(1/random);
}

int gillespieReaction(const Eigen::VectorXd& a, const double a0, const int m)
{//Calculates which reaction happens, given the rates of all reactions a(j). 
	int M=MF(m);
	double random=ud(mrand);
	double sum_aj=0;
	for (int j=0;j<M;++j){
		sum_aj += a(j);
		if(sum_aj > random*a0) {
			return j;
		}
	}
}

void gillespieUpdate(const int m, const int j, Eigen::VectorXd& x, const Eigen::MatrixXi& v)
{	//update population after reaction j happened. 
	int N=NF(m);
	for(int i=0;i<N;++i){
		x(i) += v(i,j);
	}
}

void tauLeapingk(const int m, const Eigen::VectorXd& a, const double tau, Eigen::VectorXi& k)
{//calculates k(i) the number of times reaction i happens in time tau, taken from a poisson distribution.
	int M=MF(m);
	for(int j=0;j<M;++j){
		k(j)=poisson(a(j)*tau);
	}
}

int tauLeapingUpdate(Eigen::VectorXd& xp, const Eigen::VectorXi& k,const Eigen::MatrixXi& v, const int m)
{//Updates the value of xp, the temporary population vector. returns 1 if any element of xp is negative, 0 otherwise.
int N=NF(m);
int M=MF(m);
for(int i=0;i<N;++i){
	for(int j=0;j<M;++j){
		xp(i) +=k(j)*v(i,j);
	}
	if (xp(i)<0){
		return 1;
	}
}
return 0;
}


void sim (const int m, const float xi, const float eta, const float R, const float pv0, const int exp, const int c0Exp, const int repeat, const float xinit, const float yinit, const float zinit, const float nuinit);

int main(int argc, char* argv[])
{ //main function. It takes 5 input marameters and uses them in the sim function, which does the simulation.
	if(argc!=13)
		{std::cout << "There should be 12 input parameters!\n int m, float xi, float eta, float R, float pv0, int exp, int c0Exp and int repeat.\n";}
	else{
		int m= std::atoi(argv[1]);
		float xi = std::stof(argv[2]);
		float eta = std::stof(argv[3]);
		float R = std::stof(argv[4]);
		float pv0 = std::stof(argv[5]);
		int exp = std::atoi(argv[6]);
		int c0Exp = std::stoi(argv[7]);
		int repeat = std::atoi(argv[8]);
		float xinit = std::stof(argv[9]);
		float yinit = std::stof(argv[10]);
		float zinit = std::stof(argv[11]);
		float nuinit = std::stof(argv[12]);
		sim(m,xi,eta,R,pv0,exp,c0Exp ,repeat,xinit,yinit,zinit,nuinit);
	}

}

void sim (const int m, const float xi, const float eta, const float R, const float pv0, const int exp, const int c0Exp, const int repeat, const float xinit, const float yinit, const float zinit, const float nuinit)

{//Main simulation
	clock_t t1,t2;
    	t1=clock();

	//open file for output
	std::ofstream file;
	std::string fileOut=filename(m,xi,eta,R,pv0,exp,c0Exp,repeat,xinit,yinit,zinit,nuinit);
	std::cout << "filename: \n" << fileOut << std::endl;
	file.open(fileOut);

		//Simulation parameters

	const long int nMax=pow(10,exp);//max number of tau leaping iterations
	const int gMax=500; //max number of generations simulated
	const float epsilon=0.03; //tau leaping step size
	const int maxDataPoints=4; //Save a maximum of 10^maxDataPoints data points. 
	const long int jump = jumpF(exp,maxDataPoints); //Save every jump iteration
	std::cout << "jump=" << jump << std::endl;

		//Model parameters

	const int B=170; //phage burst size
	const double c0=pow(10,c0Exp);
	//const double c0=1e9; //Number of nutrients flowing into the system
	const float g=2.4e-11; //growth rate of bacteria is g*C
	const float f=0.4*g*c0; //rate of flow in and out of the system
	const float alpha=2e-10; //rate of phage adsorption
	const float e = 1.-xi; //Spacer effectiveness
	const float r = R*g*c0; //Total rate of loss
	std::cout << "m, e, eta, r, c0 are: " << m << "\t" << e << "\t" << eta <<"\t" << r << "\t" << c0 << "\n";	

		//Initial condition
	const double nbi=round((1.0-nuinit)*xinit*c0);
	const double nvi=round(yinit*c0);
	const double nbsi=round(nuinit*xinit*c0);
	const double ci=round(zinit*c0);
	std::cout << "nbi, nbsi, nvi and ci are: " << nbi << "\t" << nbsi << "\t" << nvi << "\t" << ci << std::endl;

		//algorithm setup
 	const int N=NF(m); //#of species;
	const int M=MF(m); //# of reactions;
	Eigen::VectorXd a(M); //total rates of reactions
	Eigen::VectorXi k(M); //number of times reaction j happens during time interval
	Eigen::VectorXd xp(N); //temporary population
	Eigen::VectorXd x(N);  //populations nbo, nbi, nv, c
	Eigen::MatrixXi v(N,M); //change to population i due to reaction j
	int exit;
	int count;
	
	changeToPopulation(m,B,v);//calculates v
	initialPopulation(x,nbi,nbsi,nvi,ci,m);//x initial
	double t=0;
	long int n=0;

	saveToFile(x,file,m,n,t,0);//save initial state
	
	while(n<nMax && t*g*c0<gMax){
		//Check if extinct
		if(totalBacteria(x,m)==0){
			std::cout << "Bacteria extinction at n=" << n << "\n" ;
			break;
		}	
		if(x(m+1)==0){
			std::cout << "Phage extinction, sim stopped at n=" << n << "\n" ;
			break;
		}	
		//Total rates of reactions
		totalRates(x,a,m,g,f,c0,B,r,alpha,pv0,eta,e);//calculates a
		double a0=a.sum();

		//Calculate tau leaping time step
		double tau=tauTimeStep(v,a,x,m,epsilon); 

		if(tau<=10/a0){
		//If tau prime is smaller than 10/a0, less than 10 reactions happen during tau so algorimth is inefficient. Do gillespie for  100 iterations instead.
			for(int q=0;q<100;++q){
				totalRates(x,a,m,g,f,c0,B,r,alpha,pv0,eta,e);//calculates a
				a0=a.sum();
				double dt=gillespieTime(a0);//time step
				int j=gillespieReaction(a,a0,m);//reaction j happens
				gillespieUpdate(m,j,x,v);//updates x
				t=t+dt;				
				saveToFile(x,file,m,n+q+1,t,1.);//saves new values		
			}
			n=n+99;
		}

		else{
			//Tau leaing
			exit=0;
			count=0;

			while(exit==0 && count<10){
			//Calculate the new population, than check if all x(i) are positive. If not, redo this step with 1/2 the time step. 
				xp=x;//temporary population
				tauLeapingk(m,a,tau,k);//calculate k
				int check=tauLeapingUpdate(xp,k,v,m); //update temporary population, if any are negative, returns 1.

				if (check==0) {
				//all x(i) positive
					exit=1;
					x=xp;
					t=t+tau;

				}
				else{//redo with smaller time step
					tau=tau/2;
				}
				count += 1;
				if (count>=10){std::cout << "something wrong during tau leaping\n";}
			}

			if((n+1) % jump == 0){
				saveToFile(x,file,m,n+1,t,0);
			}

		}
	
	n=n+1;
	}
	std::cout << "Number of iterations run: " << n << "\n";
	file.close();
	t2=clock();
    	double diff=t2-t1;
   	std::cout<< "time computing: " << diff/CLOCKS_PER_SEC << " s\n";
}
