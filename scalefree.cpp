#include <iostream>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

#define step_number 200
#define N 40000
#define T  300//time
#define E 200		//samples for covering data 
#define E1 1	//samples for covering data
#define E2 1	//samples for average
// by choosing all samples equal to one, we will only have the graph of changes over time.

using namespace std;

const int	M = 10;//average degree is 4
		  

double		ka=1,
		beta=0.3,//probability of being inactive for active node
		alfa=-0.2,
		c_pay=1.,//intrinsic cost of being protected
		frac_act= 0.1,//initial fraction of act nodes
		frac1,
		frac2,		
		frac_protect=0.1,//initialfraction of protected nodes
		m0=0.2, //threshold limit  for not protected nodes
		m1,	//threshold limit for protected nodes
		sigma=1.2	,//incease factor for threshold limit
				pay[N],//pay off
				gama=30,
				alpha = 1/(gama - 1),
				p_change;


int		 	deg[N],//degree of each node
			bin[N],//its for binning of degree distribution (for plot the distribution well)
			node,
			neighbour0,
	  		G[T][E],//number of active site at time T in sample E
			G_p[T][E],
			G_np[T][E],
			G1[T][E],
			act,
			fact=0;

float 		p[N],//connection probability
	  			pk[N],//degree distribution
				g_act[T][E];

bool 		 	
			old[N],
			old1[N],
			 	active[N],
				protection[N],
				protection_0[N],
				protection_1[N],
			 	member[N],
			 	L[N];
	  		
	  		
vector < vector < int > > adj;
vector < int > protect;
vector<int>list1;
vector<int>list2;
vector < int > actv_p;//list of active sites at even time steps
vector < int > actv_np;
vector < int > actv;
inline double MyRand(){
    return (double) rand() / (RAND_MAX);
};


void Djstra (int x){//find giant connected component
	for (int i = 0; i < adj[x].size(); i++){
		if ( L[adj[x][i]]==false ){
			member[adj[x][i]] = true;
		}
	}
}


int main(){

	//________________OUTPUT FUNC___________________
	    const int Count = 7;              //Count of files
    string name = "a_";     //base pattern of file name
    ofstream outfstr[Count];      //creating array of 10 output file streams
    for(int i = 0; i < Count; i++) {   //open all file streams
         string name2 = name;
        char buffer[20] = "";
        sprintf(buffer,"%u",0);
        name2.append(buffer);
        sprintf(buffer,"%u",i);
        name2.append(buffer);
        name2.append(".txt");
        const char  * name3 = name2.c_str();  
        outfstr[i].open(name3);
    }
	//_______________END OUTPUT FUNC_________________
		
	int total_degree=0 ;
	
	int    linknumber 		= 0,
		   true_counter = 0,
		   hub				= N+1,
		   g			   = 0,
		   g2,
		   g3,
		   g_t		=0,
		   g1			=0,
		   gp				= 0,
		   gcc				= 0,
		   bind				= 0,
		   x1,
		   x2;
		   
	double constant,
		   density,
		   N_act,
		   prob,
		   connectionprob_ij,
		   rand_num1,
		   rand_num2,
		   control			= 1.2,
	       weight_sum 		= 0,
	       probability_sum  = 0,
	       ave_deg          = 0;

	long double mean_degree = 0.0000000;	       
	       
	adj.resize(N);
	
	srand(time(NULL));
	long  seed = 1.0;
	double ran2(long* idum);
	
	
	double tot_prob = 0 ;
	for (int i = 0; i < N ; i++){//asign weight w(i) for node i according 
		tot_prob += pow((i+1),-alpha);
	}
	cout << "TOTAL_PROB = " << tot_prob << endl;
	for (int i = 0; i < N ; i++){
		p[i] = pow((i+1),-alpha) / double(tot_prob);
		constant += p[i];
	}
	cout << "TOT = " << constant << endl;
		
	while (mean_degree != M ){
		select:
		// _________________the 1st node __________________________:

			rand_num1 = rand() / double(RAND_MAX);
			long double sum_of_prob1 = 0.0000000;

			for (int i = 0 ; i < N ; i++){	
				//sum_of_prob1 += ( pow(i,-alpha)/double(tot_prob) );
				sum_of_prob1 += p[i];
				if (rand_num1 <= sum_of_prob1){
					x1 = i ;
					break;
				}
			}
		// _________________the 2nd node __________________________:

			rand_num2 = rand() / double(RAND_MAX);
			long double sum_of_prob2 = 0.0000000 ;

			for (int i = 0 ; i < N ; i++){
				//sum_of_prob2 += ( pow(i,-alpha)/double(tot_prob) );
				sum_of_prob2 += p[i];
				if (rand_num2 <= sum_of_prob2){
					x2 = i;
					break;
				}
			}
			if ( x1==x2 )	goto	select;
			
			for (int i = 0; i < adj[x1].size(); i++){
				if ( adj[x1][i]==x2 ){
					goto	select;
				}
			}
			
			adj[x1].push_back(x2);
			adj[x2].push_back(x1);
			deg[x1]++;
			deg[x2]++;
			total_degree +=2;
			if (total_degree%1000==0){
				cout << " total degree = " << total_degree << endl;
			}
			mean_degree = total_degree*1.0 / N * 1.0;
	}//while
		
int x3;
		for (int i=0 ;i<N; i++){
		    x3=deg[i];
		    bin[x3]++;
		}
//==========================================================================
	for (int i=0 ;i<40; i++){
		    
		    cout << " bin = " << bin[i] << endl;
		}


	//_________________FIND GCC__________________
	for (int i = 0; i < N; i++){
		hub = (deg[i] > deg[hub]) ? i : hub;
	}
	cout << " hub = " << hub << endl;
	int nw = 0;
	Djstra ( hub );
	L[hub] = true;
	gcc++;	
	nw++;
	
	while ( true ){
		int ctrl = 0;
		for (int i = 0; i < N; i++){
			if ( member[i]==true ){
				Djstra ( i );
				member[i] = false;
				L[i] = true;// the node i with L[i]=true is amember of GCC
				ctrl++;
				gcc++;
			}
		}
		nw += 1;
		if ( ctrl==0 )   break;
	}
	for (int i = 0; i < N; i++){
		if ( L[i]==true ){
			gp++;// size of GCC
		}	
	}
	cout << " GCC = " << gp << endl;	
	//_________________END FIND GCC________________
	
	//_________________FES START___________________

//for(int sample2=0;sample2<50;sample2++){
//	sigma=0;
       for(int sample1=0; sample1<E1; sample1++){
//		alfa=alfa+1;
//	sigma=1;
//		beta=0.;
//		c_pay=0;
	for (int sample = 0; sample < E; sample++){
		for (int i = 0; i < N; i++){
			active[i] = false;
			protection_0[i]=false;
			protection_1[i]=false;
		}
		alfa=alfa+0.2;
//	beta=beta+0.01;
//	sigma=sigma+0.04;
//	c_pay=c_pay+0.1;
	for(int sample2=0;sample2<E2;sample2++){
	        m1=m0*sigma;
//		l_n=l_n+0.01;
		list1.clear();
		actv_p.clear();
		actv_np.clear();
		protect.clear();
	//___________________________________________________	
		act=0;
				
		g1=frac_protect*N;
		g=(frac_act)*N;
		while(g1>0){//make alist of active sites at the time "t0"
			 rand_num2 = rand() / double(RAND_MAX);
			 node=rand_num2*N;
			 if (protection[node]==false){
				  protect.push_back(node);//list of active sites
				  protection[node] = true;
				  g1--;
			}

		}


		
		while(g>0){//make alist of active sites at the time "t0"
			 rand_num2 = rand() / double(RAND_MAX);
			 node=rand_num2*N;
			 if (active[node]==false){
			   if (protection[node]==true ) actv_p.push_back(node);//list of active sites
			   if (protection[node]!=true ) actv_np.push_back(node);
			   active[node] = true;
			   g--;
			 }
		}

		frac2=0;



		for (int i = 0; i < N; i++){
		if(active[i]==true)  frac2++;

		}
		
		for(int i=0;i<N;i++){
				  if (protection[i]==true ) pay[i]=c_pay+(alfa*actv_p.size()*1.0/N);//list of active sites
				  if (protection[i]!=true ) pay[i]=(alfa*actv_np.size()*1.0/N);
		}


		for (int k = 0; k < T; k++){//time loop	

			g = 0;//number of active sites
			g1=0;
			g = actv_p.size()+ actv_np.size();// g is number of active site at even time
			g1 = protect.size();
   			G[k][sample] = g;
   			G1[k][sample] = g1;
   			G_p[k][sample] = actv_p.size();
   			G_np[k][sample] = actv_np.size();
   			if ( g==0 ){
   				for (int i = 0; i < T; i++){
   					G[i][sample] = 0;
   				}   			
 //  				cout << "there is not any more active site " << k << endl;
 //  				break;
   			}

		if(k==290) g2=g1;
		if(k==291) g3=g1;		
//		outfstr[1] << setprecision(8) << k << "\t" <<  g / double(gp*1.0) << endl;
//		outfstr[2] << setprecision(8) << k << "\t" <<  g1 / double(gp*1.0) << endl;

//		cout << "g_t " << g_t << endl;
			int q1 , site;
		protect.clear();
		actv_p.clear();
		actv_np.clear();
		actv.clear();

//		for (int l=0;l<N;l++){
		for (int i=0;i<N;i++){//parallel protection updating

		  if(deg[i]>0){
//	cout << "i =  "<<i<< "    k = " << k << endl;
			rand_num2 = rand() / double(RAND_MAX);
			node=rand_num2*adj[i].size();
			neighbour0=adj[i][node];
			p_change=1*1.0/(1.0+ exp((pay[neighbour0]-pay[i])/ka));
			rand_num2 = rand() / double(RAND_MAX);
			if(rand_num2<=p_change){
				if (protection[neighbour0]==true  ){
				   protection_1[i]=true;
				}
				if (protection[neighbour0]==false ) {
				   protection_1[i]=false;
				}

			}
			if(rand_num2>p_change){
				if (protection[i]==true ){
				   protection_1[i]=true;
				}
				if (protection[i]==false ) {
				   protection_1[i]=false;
				}

			}
		   }
  	
		}


		for (int i=0;i<N;i++){
			if(protection_1[i]==true){
			 protect.push_back(i);
			 protection[i]=true;
			}
			if(protection_1[i]==false){
			 protection[i]=false;
			}
		}
//		g1 = protect.size();
		




    		if ( g!=0 ){//parallel act updating
			for (int i = 0; i < N; i++){
				double N_act=0;
			      	for (int j = 0; j < adj[i].size(); j++){
    				 int neighbour0;
    				 neighbour0 = adj[i][j];
    				 if ( active[neighbour0]==true&& old[neighbour0]==false){
				     if(protection[neighbour0]==false ) N_act++;
				     if(protection[neighbour0]==true) N_act=N_act+0.5;
				 }
    						
    				} 
				rand_num2 = rand()/double(RAND_MAX);
				if (protection[i]==false &&(N_act*1.0/(deg[i]*1.0))>=m0 && active[i]==false){ 
				   active[i]=true;
				   old[i]=true;
				}
				if (protection[i]==true &&(N_act*1.0/(deg[i]*1.0))>=m1 && active[i]==false){ 
				   active[i]=true;
				   old[i]=true;
				}
			}

		}
		for(int i=0;i<N;i++){
				rand_num2 = rand()/double(RAND_MAX);
				if ( active[i]==true&&  rand_num2<=beta && old[i]==false){ 
				   active[i]=false;
				}
				old[i]=false;
		}//end parallel act updating


		for (int i=0;i<N;i++){
			if(active[i]==true && protection[i]==true ) actv_p.push_back(i);
			if(active[i]==true && protection[i]==false ) actv_np.push_back(i);    
			   
		}


		for (int i=0;i<N;i++){
			if(protection[i]==true)    pay[i]=c_pay+alfa*(actv_p.size()*1.0/N);
			else pay[i]=alfa*(actv_np.size()*1.0/N);
		}

		for (int i=0;i<N;i++){
			if(active[i]==true) actv.push_back(i);   
		}


//		outfstr[1] << setprecision(8) << k << "\t" <<  g / double(gp*1.0) << endl;
//		outfstr[2] << setprecision(8) << k << "\t" <<  g1 / double(gp*1.0) << endl;


	}//for k

	}//E2

	if((g*1.0/N)>0.2) fact++;
		
		int g_t_p=0;

		  for(int i=0;i<400;i++){//time average
		     g_t_p=g_t_p+G_p[T-i][sample];
		  }


		int g_t_np=0;

		  for(int i=0;i<400;i++){//time average
		     g_t_np=g_t_np+G_np[T-i][sample];
		  }


		

	
		int g_t=0;

		  for(int i=0;i<400;i++){//time average
		     g_t=g_t+G[T-i-10][sample];
		  }
		
		int g1_t=0;

		  for(int i=0;i<400;i++){
		     g1_t=g1_t+G1[T-i-10][sample];
		  }

		outfstr[1] << setprecision(8) << alfa << "\t" << G[290][sample]*1.0 / double(gp*1.0) << endl;
		outfstr[2] << setprecision(8) << alfa << "\t" << G1[290][sample]*1.0 / double(gp*1.0)<< endl;
		outfstr[2] << setprecision(8) << alfa << "\t" << G1[291][sample]*1.0 / double(gp*1.0) << endl;
		outfstr[1] << setprecision(8) << alfa << "\t" << G[291][sample]*1.0 / double(gp*1.0) << endl;
 //  		actv_site.clear();
//		outfstr[1] << setprecision(8) << "\t" <<  g_t*1.0 / double(gp*1.0*400) << "\t";
//		outfstr[2] << setprecision(8) << "\t" <<  g1_t*1.0 / double(gp*1.0*400) <<"\t";
//		outfstr[5] << setprecision(8) << "\t" <<  g_t_p*1.0 / double(gp*1.0*400) <<"\t";
//		outfstr[6] << setprecision(8) << "\t" <<  g_t_np*1.0 / double(gp*1.0*400) <<"\t";
   }//for E
		outfstr[1] << setprecision(8) << endl;
		outfstr[2] << setprecision(8) << endl;
//		outfstr[5] << setprecision(8) << endl;
//		outfstr[6] << setprecision(8) << endl;
   		cout << " sample => " << sample1 <<"  g =  "<<g_t<< endl;
  }//E1
/*	      for(int i=0;i<5;i++){
		outfstr[3] << setprecision(8) << endl;
		outfstr[4] << setprecision(8) << endl;
		outfstr[5] << setprecision(8) << endl;
		outfstr[6] << setprecision(8) << endl;
	      }
   		cout << " sample2 => " << sample2 <<"  g =  "<<g<< endl;
//}//E2	
*/

  
	float lastSum = 0;
  

	cout << "final sum  = " << fact << endl;
	outfstr[0] << "Gcc = " << gp << endl;

	return 0;		
}
//===============ran2===============
double ran2(long *idum) {
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

	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;





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
}
