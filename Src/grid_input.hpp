#ifndef GPLUTO_GRID_INPUT_HPP
#define GPLUTO_GRID_INPUT_HPP
#include "grid_spec.hpp"
#include <algorithm>
#include <cmath>
#include <climits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace GridInput {
inline int Integer(const std::string &s) {
  size_t n; long v=std::stol(s,&n);
  if(n!=s.size() || v<1 || v>INT_MAX) throw std::runtime_error("expected a positive integer: "+s);
  return (int)v;
}
inline double Real(const std::string &s) {
  size_t n; double v=std::stod(s,&n);
  if(n!=s.size() || !std::isfinite(v)) throw std::runtime_error("expected a finite coordinate: "+s);
  return v;
}
inline void Parse(const std::string &line, GridSpec &spec,
                  double &left,int &total,double &right,bool active=true) {
  spec=GridSpec{};
  std::istringstream input(line.substr(0,line.find('#')));
  std::string label,s; input>>label;
  std::vector<std::string> t;
  while(input>>s) t.push_back(s);
  if(t.size()==3) {
    left=Real(t[0]); total=Integer(t[1]); right=Real(t[2]);
    if(!(right>left) && (active || right!=left)) throw std::runtime_error("grid endpoints must increase");
    return;
  }
  if(t.empty()) throw std::runtime_error("empty grid specification");
  spec.count=Integer(t[0]);
  if(spec.count>GRID_MAX_PATCHES || t.size()!=size_t(3*spec.count+2))
    throw std::runtime_error("expected Npatch x0 N0 type0 x1 ... (at most 16 patches)");
  spec.edge[0]=Real(t[1]); total=0;
  for(int p=0;p<spec.count;++p) {
    int n=spec.cells[p]=Integer(t[2+3*p]);
    if(total>INT_MAX-n) throw std::runtime_error("grid cell count overflow");
    total+=n;
    s=t[3+3*p];
    if(s=="u" || s=="uniform") spec.type[p]='u';
    else if(s=="s" || s=="stretched") spec.type[p]='s';
    else if(s=="l+") spec.type[p]='+';
    else if(s=="l-") spec.type[p]='-';
    else throw std::runtime_error("unknown grid type: "+s);
    spec.edge[p+1]=Real(t[4+3*p]);
    if(!(spec.edge[p+1]>spec.edge[p]) &&
       (active || spec.count!=1 || spec.type[p]!='u' || spec.edge[p+1]!=spec.edge[p])) throw std::runtime_error("patch endpoints must increase");
    if((spec.type[p]=='+' || spec.type[p]=='-') && spec.edge[p]==0.)
      throw std::runtime_error("logarithmic patch left endpoint must be nonzero");
  }
  left=spec.edge[0]; right=spec.edge[spec.count];
}
// Stable sum r+r^2+...+r^n; monotone on r>0, including contracting patches.
inline double Sum(double r,int n) {
  if(std::abs(r-1.)<1.e-12) return (double)n;
  return r*std::expm1(n*std::log(r))/(r-1.);
}
inline std::vector<double> Edges(const GridSpec &s,int factor=1) {
  if(s.count<1 || factor<1) throw std::runtime_error("invalid grid specification/coarsening");
  int n[GRID_MAX_PATCHES],beg[GRID_MAX_PATCHES];
  bool done[GRID_MAX_PATCHES]={}; int total=0,remaining=s.count;
  for(int p=0;p<s.count;++p) {
    if(s.cells[p]%factor) throw std::runtime_error("coarsening must divide every patch cell count");
    beg[p]=total; n[p]=s.cells[p]/factor; total+=n[p];
  }
  std::vector<double> x(total+1),width(total);
  for(int p=0;p<s.count;++p) {
    if(s.type[p]=='s') continue;
    double l=s.edge[p],r=s.edge[p+1],a=std::abs(l);
    double logratio=s.type[p]=='u'?0.:std::log1p((r-l)/a);
    for(int j=0;j<=n[p];++j) {
      double q=double(j)/n[p];
      x[beg[p]+j]=s.type[p]=='u'?l+(r-l)*q:
        s.type[p]=='+'?l+a*std::expm1(q*logratio):
                      r-a*std::expm1((1.-q)*logratio);
    }
    x[beg[p]]=l; x[beg[p]+n[p]]=r;
    for(int j=0;j<n[p];++j) width[beg[p]+j]=x[beg[p]+j+1]-x[beg[p]+j];
    done[p]=true; --remaining;
  }
  while(remaining) {
    bool progress=false;
    for(int p=0;p<s.count;++p) {
      if(done[p]) continue;
      int neighbor=-1;
      if(p+1<s.count && done[p+1]) neighbor=p+1;
      if(p>0 && done[p-1]) neighbor=p-1; // PLUTO prefers the left neighbor.
      if(neighbor<0) continue;
      double w=neighbor<p?width[beg[p]-1]:width[beg[p]+n[p]];
      double target=(s.edge[p+1]-s.edge[p])/w;
      double lo=0.,hi=2.;
      while(Sum(hi,n[p])<target) hi*=2.;
      for(int it=0;it<120;++it) {
        double mid=.5*(lo+hi);
        if(Sum(mid,n[p])<target) lo=mid; else hi=mid;
      }
      double ratio=.5*(lo+hi);
      x[beg[p]]=s.edge[p]; x[beg[p]+n[p]]=s.edge[p+1];
      for(int j=0;j<n[p];++j) {
        w*=ratio;
        int i=neighbor<p?beg[p]+j:beg[p]+n[p]-1-j;
        width[i]=w;
        if(neighbor<p) x[i+1]=x[i]+w; else x[i]=x[i+1]-w;
      }
      x[beg[p]]=s.edge[p]; x[beg[p]+n[p]]=s.edge[p+1];
      done[p]=true; --remaining; progress=true;
      break; // Match PLUTO's scan order for chained stretched patches.
    }
    if(!progress) throw std::runtime_error("stretched grid needs a uniform or logarithmic anchor patch");
  }
  for(int i=0;i<total;++i)
    if(!std::isfinite(x[i]) || !std::isfinite(x[i+1]) || !(x[i+1]>x[i]))
      throw std::runtime_error("grid contains nonfinite or non-increasing edges");
  return x;
}
}
#endif
