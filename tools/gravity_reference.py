#!/usr/bin/env python3
"""Independent spherical-harmonic acceleration fixtures from DE440 comments.

Decimal Legendre-polynomial potentials and five-point derivatives; no Spp code.
"""
import argparse
import csv
from decimal import Decimal as D, getcontext
from fractions import Fraction as F
import hashlib
from pathlib import Path
import re
import struct
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel',type=Path)
    args=parser.parse_args()
    raw=args.kernel.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!='c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2':
        raise ValueError('review kernel revision before regenerating force constants')
    record=struct.unpack_from('<i',raw,76)[0]
    comments=b''.join(raw[i:i+1000] for i in range(1024,(record-1)*1024,1024)).replace(b'\0',b'\n').decode('ascii')
    def parameter(name):
        match=re.search(r'(?m)^'+name+r'\s+([-+\d.DEded]+)',comments)
        if not match:raise ValueError('missing constant '+name)
        return D(match[1].replace('D','E').replace('d','e'))
    getcontext().prec=70
    # Rational ordinary Legendre polynomials. Their mth derivatives generate
    # the unnormalized associated functions without a phase convention sign.
    polynomials=[[F(1)],[F(0),F(1)]]
    for n in range(2,7):
        p=[F(0)]*(n+1)
        for i,a in enumerate(polynomials[-1]):p[i+1]+=F(2*n-1,n)*a
        for i,a in enumerate(polynomials[-2]):p[i]-=F(n-1,n)*a
        polynomials.append(p)
    derivatives={}
    for n,p in enumerate(polynomials):
        for m in range(n+1):
            derivatives[n,m]=[D(v.numerator)/D(v.denominator) for v in p]
            p=[(i+1)*v for i,v in enumerate(p[1:])]
    mu=D('403503.2356254802')*D(86400)**2
    rows=[]
    for body,max_degree,radius in [('earth',5,parameter('RE')),('moon',6,parameter('AM'))]:
        c,s={},{}
        for n in range(2,max_degree+1):
            c[n,0]=-parameter(f'J{n}'+('E' if body=='earth' else 'M'))
            if body=='moon':
                for m in range(1,n+1):
                    c[n,m]=parameter(f'C{n}{m}M');s[n,m]=parameter(f'S{n}{m}M')
        def potential(point,degree):
            x,y,z=point;r=sum(v*v for v in point).sqrt();sin=z/r
            real,imag=D(1),D(0);angles=[]
            for m in range(degree+1):
                angles.append((real,imag));real,imag=real*x-imag*y,real*y+imag*x
            total=D(0)
            for n in range(2,degree+1):
                for m in range(n+1):
                    p=sum(a*(sin**i if i else D(1)) for i,a in enumerate(derivatives[n,m]))
                    real,imag=angles[m]
                    total+=(radius/r)**n*p*(c.get((n,m),D(0))*real+s.get((n,m),D(0))*imag)/r**m
            return mu/r*total
        points=[('384400','0','0'),('0','-370000','0'),('0','0','400000'),('0','0','-350000'),
                ('-332754.267845','179885.015158','-2950.415724'),('200000','-300000','70000'),
                (str(radius*D('1.2')),str(radius*D('0.2')),str(radius*D('0.3')))]
        for degree in (2,max_degree):
            for values in points:
                point=list(map(D,values));acceleration=[];h=D('0.001')
                for axis in range(3):
                    f=[]
                    for offset in (-2,-1,1,2):
                        q=point.copy();q[axis]+=offset*h;f.append(potential(q,degree))
                    acceleration.append((f[0]-8*f[1]+8*f[2]-f[3])/(12*h))
                rows.append([body,degree,*values,str(radius),str(mu),*[format(a,'.25e') for a in acceleration]])
    def write(name,header,values):
        target=Path(__file__).resolve().parent/'reference'/name;tmp=None
        try:
            with tempfile.NamedTemporaryFile('w',dir=target.parent,delete=False,newline='') as out:
                tmp=Path(out.name);out.write('# Independent Decimal potential gradients; DE440s SHA256 '+digest+'\n'+header)
                csv.writer(out,lineterminator='\n').writerows(values)
            tmp.replace(target)
        finally:
            if tmp is not None:tmp.unlink(missing_ok=True)
        print(f'{target}: {len(values)} states')
    write('gravity_de440.csv','# body,degree,x,y,z,R_km,mu_km3_day2,ax,ay,az_km_day2\n',rows)
    tide_rows=[];radius=parameter('RE');mu=D('4902.8001184575496')*D(86400)**2
    # Moon GM is the published gm_de440.tpc value, converted to km^3/day^2.
    # The potential below has
    # independent scalar m=0,1,2 expansions about the Z pole.
    for order in range(3):
        for values in [('384400','0','0'),('200000','-300000','70000'),('0','0','380000')]:
            for source in [('210000','180000','-230000'),('1.4e8','-2e7','5e7')]:
                point=list(map(D,values));q=list(map(D,source));love=parameter(f'K2E{order}')
                def potential(p):
                    x,y,z=p;qx,qy,qz=q;r2=x*x+y*y+z*z;q2=qx*qx+qy*qy+qz*qz
                    if order==0:shape=(qz*qz-q2/3)*(D('1.5')*z*z-r2/2)
                    elif order==1:shape=2*qz*z*(x*qx+y*qy)
                    else:shape=(x*qx+y*qy)**2-(x*x+y*y)*(qx*qx+qy*qy)/2
                    return D('1.5')*love*mu*radius**5*shape/(q2*q2*q2.sqrt()*r2*r2*r2.sqrt())
                acceleration=[];h=D('0.001')
                for axis in range(3):
                    f=[]
                    for offset in (-2,-1,1,2):
                        p=point.copy();p[axis]+=offset*h;f.append(potential(p))
                    acceleration.append((f[0]-8*f[1]+8*f[2]-f[3])/(12*h))
                tide_rows.append([order,*values,*source,str(love),str(radius),str(mu),*[format(a,'.25e') for a in acceleration]])
    write('earth_tides.csv','# order,rx,ry,rz,qx,qy,qz,love,R_km,mu_km3_day2,ax,ay,az_km_day2\n',tide_rows)


if __name__=='__main__':main()
