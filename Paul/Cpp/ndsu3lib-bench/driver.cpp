// ndsu3lib benchmark driver: (7,7) x (7,7) -> each of the 113 irreps of GT [14,7,0]^2.
// Usage: bench_driver irreps.txt time [reps] | dump out.txt | mult
// GT mu=[m1,m2,m3] of the coupled irrep <-> ndsu3lib (lambda,mu)_nd = (m2-m3, m1-m2)  (nd HW state has eps=2*lam+mu, Lambda=mu/2,
// which matches the GT HW state Lambda=(m1-m2)/2). Sites: (7,7) both (self-conjugate, no ambiguity).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "ndsu3lib.h"
using namespace std;
struct Ir { int m1,m2,m3,m; ndsu3lib::su3irrep nd; int rhomax; };
static vector<Ir> load(const char* fn){
  vector<Ir> v; FILE* f=fopen(fn,"r"); char line[256];
  while(fgets(line,256,f)){ if(line[0]=='#')continue; Ir r; if(sscanf(line,"%d %d %d %d",&r.m1,&r.m2,&r.m3,&r.m)==4){
      r.nd.lambda=r.m2-r.m3; r.nd.mu=r.m1-r.m2; v.push_back(r);} }
  fclose(f); return v; }
static const ndsu3lib::su3irrep S={7,7};
static vector<int> P1,P2,Q2; static vector<double> W;
static int call(const Ir& r,int eps3,int L32,int& numb){
  int dimpq=(max(S.lambda,S.mu)+1)*(S.lambda+1)*(S.mu+1);
  P1.resize(dimpq);P2.resize(dimpq);Q2.resize(dimpq);W.resize((size_t)dimpq*r.rhomax);
  ndsu3lib::calculate_coupling_canonical(S,S,r.nd,eps3,L32,dimpq,dimpq*r.rhomax,r.rhomax,numb,W.data(),P1.data(),P2.data(),Q2.data());
  return numb*r.rhomax; }
// all (eps3,2*Lambda3) states of irrep3: GT pairs (m12,m22), l1>=m12>=l2>=m22>=l3
template<class F> static void states(const Ir& r,F f){
  for(int a=r.m2;a<=r.m1;a++)for(int b=r.m3;b<=r.m2;b++){ int eps=3*(a+b)-2*(r.m1+r.m2+r.m3); f(eps,a-b);} }
static double now(){ return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc,char**argv){
  vector<Ir> irr=load(argv[1]); string mode=argv[2];
  ndsu3lib::initialize_ndsu3lib(false,false,200);
  int bad=0; for(auto&r:irr){ r.rhomax=ndsu3lib::outer_multiplicity(S,S,r.nd); if(r.rhomax!=r.m){bad++;printf("MULT MISMATCH [%d %d %d] nd(%d,%d) rho=%d m=%d\n",r.m1,r.m2,r.m3,r.nd.lambda,r.nd.mu,r.rhomax,r.m);} }
  if(mode=="mult"){ int s=0; for(auto&r:irr)s+=r.rhomax; printf("irreps %zu, sum rhomax %d, mismatches %d\n",irr.size(),s,bad); return 0; }
  if(mode=="dump"){ FILE*o=fopen(argv[3],"w");
    for(auto&r:irr){ int eps3=r.m1+r.m2-2*r.m3, numb; call(r,eps3,r.nd.mu,numb);
      fprintf(o,"IRREP %d %d %d %d %d\n",r.m1,r.m2,r.m3,r.rhomax,numb);
      for(int i=0;i<numb;i++){ if(P2[i]!=0||Q2[i]!=0)continue; // site 2 at its highest weight (p2=q2=0)
        int p1=P1[i]; int q1=(2*(S.lambda+S.lambda)+S.mu+S.mu-eps3)/3-p1; int e1=eps3-(2*S.lambda+S.mu); int L1=S.mu+p1-q1;
        fprintf(o,"C %d %d",e1,L1); for(int rho=0;rho<r.rhomax;rho++)fprintf(o," %.17g",W[i+numb*rho]); fprintf(o,"\n"); } }
    fclose(o); return 0; }
  int reps=argc>3?atoi(argv[3]):5;
  // (a) full set (every eps3,Lambda3 of every irrep); (b) HW row only (eps3max, Lambda3=mu_nd/2)
  for(int rep=0;rep<reps;rep++){
    long ncf=0,ncalls=0; double t0=now();
    for(auto&r:irr) states(r,[&](int e,int L){int nb; ncf+=call(r,e,L,nb); ncalls++;});
    double ta=now()-t0;
    long nhw=0; t0=now();
    for(auto&r:irr){int nb; nhw+=call(r,r.m1+r.m2-2*r.m3,r.nd.mu,nb);}
    double tb=now()-t0;
    printf("rep %d  full: %.4f s (%ld calls, %ld coefficients)   HWrow: %.5f s (113 calls, %ld coefficients)\n",rep,ta,ncalls,ncf,tb,nhw);
  }
  // single-irrep timing for the largest one
  return bad?1:0; }
