#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>

#include "root_refine.hpp"

namespace spp::detail {

struct ScalarRoot {
    double x;
    bool rising;
    bool contact;
};

// Isolate roots on [0,width] using a supplied bound on |f''|. The deviation
// from the endpoint chord is at most curvature*h*h/8. Unlike a sign-only
// scan, subdivision can therefore find two crossings with same-sign endpoints
// and stationary contacts. The bound is a caller/model contract, checked at
// each midpoint; it is not inferred from a sparse sampling of an arbitrary f.
template<class F>
std::vector<ScalarRoot> bounded_roots(F function,double width,double tolerance,
                                     double contact_tolerance,double curvature,
                                     double roundoff=1e-12) {
    if(!std::isfinite(width) || width<=0 || !std::isfinite(tolerance) || tolerance<=0
        || !std::isfinite(contact_tolerance) || contact_tolerance<0
        || !std::isfinite(curvature) || curvature<0 || !std::isfinite(roundoff) || roundoff<0)
        throw std::invalid_argument("event search: invalid interval, tolerance or curvature");
    int evaluations=0;
    auto f=[&](double x) {
        if(++evaluations>65536)
            throw std::runtime_error("event search unresolved: continuous contact or incompatible model");
        const double y=function(x);
        if(!std::isfinite(y))throw std::runtime_error("event search returned a non-finite value");
        return y;
    };
    auto refine=[&](double a,double fa,double b,double fb) {
        if(fa==0)return a;
        if(fb==0)return b;
        return fa>0?refine_root(f,a,fa,b,fb,tolerance,100)
                   :refine_root(f,b,fb,a,fa,tolerance,100);
    };
    struct Contact {double a,b,best,error;};
    std::vector<ScalarRoot> roots;
    std::vector<Contact> contacts;
    auto walk=[&](auto&& self,double a,double fa,double b,double fb,int depth)->void {
        const double h=b-a;
        const double error=curvature*h*h/8+roundoff;
        if(std::min(fa,fb)>error+contact_tolerance
            || std::max(fa,fb)<-error-contact_tolerance)return;
        const double m=a+0.5*h,fm=f(m);
        if(std::abs(fm-0.5*(fa+fb))>error+roundoff)
            throw std::runtime_error(std::format("event search exceeds curvature bound or is discontinuous: [{:.17g}, {:.17g}], chord deviation {:.9g}, bound {:.9g}",
                a,b,std::abs(fm-0.5*(fa+fb)),error+roundoff));
        const bool crossing=(fa<0 && fb>0)||(fa>0 && fb<0);
        const bool monotonic=std::abs(fb-fa)>curvature*h*h/2+2*roundoff;
        if(monotonic) {
            if(crossing || fa==0 || fb==0)roots.push_back({refine(a,fa,b,fb),fb>fa,false});
            return;
        }
        if(h<=tolerance || depth>=48) {
            if(crossing)roots.push_back({refine(a,fa,b,fb),fb>fa,false});
            else {
                double x=m,y=std::abs(fm);
                if(std::abs(fa)<y){x=a;y=std::abs(fa);}
                if(std::abs(fb)<y){x=b;y=std::abs(fb);}
                if(y<=contact_tolerance+roundoff)contacts.push_back({a,b,x,y});
            }
            return;
        }
        self(self,a,fa,m,fm,depth+1);self(self,m,fm,b,fb,depth+1);
    };
    walk(walk,0,f(0),width,f(width),0);

    // Merge the small unresolved neighborhoods around one stationary contact.
    // Real crossings found within a neighborhood take precedence; retain both
    // rise and set if a short above-horizon interval really exists.
    for(std::size_t i=0;i<contacts.size();) {
        Contact c=contacts[i++];
        while(i<contacts.size() && contacts[i].a<=c.b+2*tolerance) {
            if(contacts[i].error<c.error){c.best=contacts[i].best;c.error=contacts[i].error;}
            c.b=contacts[i++].b;
        }
        bool has_crossing=false;
        for(const auto& r:roots)if(r.x>=c.a-2*tolerance && r.x<=c.b+2*tolerance){has_crossing=true;break;}
        if(has_crossing)continue;
        const double pad=std::max(10*tolerance,c.b-c.a);
        double a=std::max(0.0,c.a-pad),b=std::min(width,c.b+pad);
        const double fa=std::abs(f(a)),fb=std::abs(f(b));
        if(c.a<=tolerance && std::abs(f(0))<=contact_tolerance) {
            const double d=std::min(width/2,std::max(1.0/86400,10*tolerance));
            const double derivative=(-3*f(0)+4*f(d)-f(2*d))/(2*d);
            if(std::abs(derivative)<=2*curvature*tolerance+4*roundoff/d) {
                roots.push_back({0,false,true});continue;
            }
        }
        if(c.b>=width-tolerance && std::abs(f(width))<=contact_tolerance) {
            const double d=std::min(width/2,std::max(1.0/86400,10*tolerance));
            const double derivative=(3*f(width)-4*f(width-d)+f(width-2*d))/(2*d);
            if(std::abs(derivative)<=2*curvature*tolerance+4*roundoff/d) {
                roots.push_back({width,false,true});continue;
            }
        }
        if(c.best<=a || c.best>=b || !(c.error<fa && c.error<fb))continue;
        // Minimize the signed function towards zero, rather than |f|. The
        // latter has a minimum at every crossing and can invent an additional
        // contact immediately beside a real rise/set pair.
        const double sign=f(a)<0?-1.0:1.0;
        auto objective=[&](double t){return sign*f(t);};
        constexpr double ratio=0.6180339887498948482;
        double x=b-ratio*(b-a),y=a+ratio*(b-a),fx=objective(x),fy=objective(y);
        for(int j=0;j<100 && b-a>tolerance;++j) {
            if(fx<fy){b=y;y=x;fy=fx;x=b-ratio*(b-a);fx=objective(x);}
            else {a=x;x=y;fx=fy;y=a+ratio*(b-a);fy=objective(y);}
        }
        const double best=fx<fy?x:y;
        // A contact needs a stationary interior extremum. A monotonic flank
        // within the angular tolerance is not an additional event.
        const double guard=std::max(tolerance,std::sqrt(8*roundoff/std::max(curvature,1e-12)));
        const double value=objective(best);
        if(best>guard && best+guard<width && value>=-roundoff
            && std::abs(value)<=contact_tolerance+roundoff
            && objective(best-guard)>=value && objective(best+guard)>=value)
            roots.push_back({best,false,true});
    }
    const double guard=std::max(tolerance,std::sqrt(8*roundoff/std::max(curvature,1e-12)));
    for(auto& root:roots)if(!root.contact && root.x>guard && root.x+guard<width) {
        const double left=f(root.x-guard),right=f(root.x+guard);
        if((left>0 && right>0)||(left<0 && right<0))root.contact=true;
    }
    std::sort(roots.begin(),roots.end(),[](const auto& a,const auto& b){return a.x<b.x;});
    std::vector<ScalarRoot> result;
    for(const auto& root:roots) {
        if(!result.empty() && root.contact && result.back().contact && root.x-result.back().x<=2*guard) {
            result.back().x=0.5*(result.back().x+root.x);continue;
        }
        if(!result.empty() && root.x-result.back().x<=tolerance
            && root.contact==result.back().contact && root.rising==result.back().rising)continue;
        result.push_back(root);
    }
    return result;
}

}  // namespace spp::detail
