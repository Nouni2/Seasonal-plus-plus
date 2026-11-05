#!/usr/bin/env python3
"""Generate independent IAU-model and JPL lunar-rotation benchmark fixtures.

No Spp code is loaded or executed. Normal builds use the committed CSVs.
"""
import argparse
import csv
from decimal import Decimal as D
import hashlib
import json
import math
from pathlib import Path
import re
import tempfile
import urllib.parse
import urllib.request


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    dest=Path(__file__).resolve().parent/'reference'

    def fetch(name,url):
        path=args.cache/(name+'-'+hashlib.sha256(url.encode()).hexdigest()[:16])
        if not path.exists():
            with urllib.request.urlopen(url,timeout=45) as response:
                content=response.read()
            path.write_bytes(content)
        return path.read_bytes()

    def write(name,header,rows):
        temporary=None
        try:
            with tempfile.NamedTemporaryFile('w',newline='',dir=dest,suffix='.tmp',delete=False) as file:
                temporary=Path(file.name)
                file.write(header)
                csv.writer(file,lineterminator='\n').writerows(rows)
            temporary.replace(dest/name)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    pck=fetch('pck00011.tpc','https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc')
    text=pck.decode('ascii')
    def constants(name):
        group=re.search(r'^\s*'+name+r'\s*=\s*\((.*?)\)',text,re.M|re.S)
        if not group:
            raise ValueError('missing PCK field '+name)
        return [D(x.replace('D','E').replace('d','e')) for x in group[1].split()]
    a,b,w=[constants('BODY301_'+field) for field in ('POLE_RA','POLE_DEC','PM')]
    aa,bb,ww=[constants('BODY301_NUT_PREC_'+field) for field in ('RA','DEC','PM')]
    angles=constants('BODY3_NUT_PREC_ANGLES')
    def model(jd):
        days=D(str(jd))-D('2451545');T=days/D('36525')
        alpha=float(a[0]+a[1]*T+a[2]*T*T)
        delta=float(b[0]+b[1]*T+b[2]*T*T)
        meridian=float((w[0]+w[1]*days+w[2]*days*days)%360)
        for i in range(13):
            angle=math.radians(float((angles[2*i]+angles[2*i+1]*T)%360))
            alpha+=float(aa[i])*math.sin(angle)
            delta+=float(bb[i])*math.cos(angle)
            meridian+=float(ww[i])*math.sin(angle)
        alpha%=360;meridian%=360
        r,d,m=map(math.radians,(alpha,delta,meridian))
        pole=[math.cos(d)*math.cos(r),math.cos(d)*math.sin(r),math.sin(d)]
        node=[-math.sin(r),math.cos(r),0]
        north=[-math.sin(d)*math.cos(r),-math.sin(d)*math.sin(r),math.cos(d)]
        first=[math.cos(m)*u+math.sin(m)*v for u,v in zip(node,north)]
        second=[-math.sin(m)*u+math.cos(m)*v for u,v in zip(node,north)]
        return [jd,alpha,delta,meridian,*first,*second,*pole]
    epochs=[1319270.0,2305447.5,2415020.5,2451545.0,2460409.25,2461348.0,2488069.5,2853320.0]
    write('moon_iau_rotation.csv',
          '# Independently parsed NAIF pck00011.tpc; decimal phase arithmetic; TDB.\n'
          '# PCK SHA256 '+hashlib.sha256(pck).hexdigest()+'\n'
          '# jd_tdb,pole_ra,pole_dec,prime_meridian,r00,r01,r02,r10,r11,r12,r20,r21,r22\n',
          [model(jd) for jd in epochs])

    # High-precision lunar ME orientation: recover axes from independent
    # astrometric line of sight, sub-observer longitude/latitude and pole.
    times=[2415020.5,2451545.0,2460409.25,2461348.0,2461369.0,2469807.5,2488069.5]
    params=dict(format='json',COMMAND="'301'",OBJ_DATA="'NO'",EPHEM_TYPE="'OBSERVER'",
                CENTER="'500@399'",TIME_TYPE="'TT'",TLIST_TYPE="'JD'",
                TLIST="'"+' '.join(map(str,times))+"'",QUANTITIES="'1,14,21,32'",
                REF_SYSTEM="'ICRF'",ANG_FORMAT="'DEG'",APPARENT="'AIRLESS'",
                EXTRA_PREC="'YES'",CSV_FORMAT="'YES'",CAL_FORMAT="'BOTH'")
    value=json.loads(fetch('rotation-de441.json','https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(params)))
    if value.get('error'):
        raise ValueError(value['error'])
    result=value['result']
    if 'DE441' not in result or 'MOON_ME' not in result or 'Date__(TT)' not in result:
        raise ValueError('review reference model/time/frame change')
    rows=[]
    for row in csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()):
        # date, JD, solar/lunar flags, ICRF RA/Dec, sub-observer lon/lat,
        # one-way light time in minutes, ICRF north-pole RA/Dec.
        rows.append([float(row[i]) for i in (1,8,4,5,6,7,9,10)])
    if len(rows)!=len(times):
        raise ValueError('reference epoch count changed')
    write('moon_rotation_horizons.csv',
          '# JPL DE441 Moon ME: TT reception epoch and one-way light time.\n'
          '# jd_tt,light_minutes,icrf_ra,icrf_dec,sub_lon,sub_lat,pole_ra,pole_dec\n',rows)
    surface_rows=[]
    for name,lat,lon,elev,epochs in [('geocenter',0,0,0,times),
            ('paris',48.8566,2.3522,35,[2460409.25,2461348.0,2461369.0]),
            ('sydney',-33.8688,151.2093,58,[2460409.25,2461348.0,2461369.0])]:
        request=dict(params,QUANTITIES="'14,15,17,49'",TLIST="'"+' '.join(map(str,epochs))+"'")
        if name!='geocenter':
            request.update(CENTER="'coord@399'",COORD_TYPE="'GEODETIC'",SITE_COORD=f"'{lon},{lat},{elev/1000}'")
        data=json.loads(fetch('surface-'+name+'.json','https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(request)))
        if data.get('error'):
            raise ValueError(data['error'])
        response=data['result']
        if 'DE441' not in response or 'MOON_ME' not in response or 'Date__(TT)' not in response:
            raise ValueError('surface reference model/time/frame changed')
        values=list(csv.reader(response.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))
        if len(values)!=len(epochs):
            raise ValueError('surface reference count changed')
        for row in values:
            number=lambda i: float(row[i]) if row[i].strip()!='n.a.' else 'nan'
            surface_rows.append([name,number(1),lat,lon,elev,number(10),
                                 number(4),number(5),number(6),number(7),number(8)])
    write('moon_surface_horizons.csv',
          '# JPL DE441 MOON_ME, TT reception, geocentric or geodetic observer.\n'
          '# source,jd_tt,lat,lon,elev_m,dut1,sub_lon,sub_lat,solar_lon,solar_lat,pole_pa\n',surface_rows)
    erfa=fetch('erfa-tests.c','https://raw.githubusercontent.com/liberfa/erfa/master/src/t_erfa_c.c')
    block=re.search(r'static void t_bp06\(.*?(?=static void)',erfa.decode(),re.S)
    if not block:
        raise ValueError('ERFA bp06 reference block missing')
    bias=[]
    for i in range(3):
        row=[]
        for j in range(3):
            value=re.search(r'vvd\(rb\['+str(i)+r'\]\['+str(j)+r'\],\s*([-+\d.eE]+)',block[0])
            if not value:
                raise ValueError('ERFA frame-bias coefficient missing')
            row.append(float(value[1]))
        bias.append(row)
    write('frame_bias_erfa.csv',
          '# ICRF to mean equatorial J2000 bias matrix, ERFA t_bp06 independent reference.\n'
          '# Source SHA256 '+hashlib.sha256(erfa).hexdigest()+'\n',bias)
    print('Wrote independent IAU and DE441 rotation benchmark references')


if __name__=='__main__':
    main()
