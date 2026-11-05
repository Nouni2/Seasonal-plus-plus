#!/usr/bin/env python3
"""Generate SPK-reader reference states with independent CSPICE (SpiceyPy).

Supply an existing DE440s kernel. This script never loads or executes Spp.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import tempfile
import spiceypy as spice


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel',type=Path)
    args=parser.parse_args()
    digest=hashlib.sha256(args.kernel.read_bytes()).hexdigest()
    if digest!='c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2':
        raise ValueError('review reference change: expected published DE440s kernel')
    spice.furnsh(str(args.kernel))
    cover=spice.spkcov(str(args.kernel),301)
    first,last=spice.wnfetd(cover,0)
    # Exact ends, internal record boundaries and non-binary sub-day instants.
    epochs=[first,first+1,first+345600-1e-4,first+345600,0.0,
            764445600.123456,846244800,846288000.321,847497600,
            last-345600,last-1,last]
    pairs=[(301,399),(399,301),(301,3),(399,3),(3,0),(10,399),
           (1,399),(2,399),(4,399),(5,399),(6,399),(7,399),(8,399),(9,399),(399,0)]
    rows=[]
    for et in epochs:
        for target,center in pairs:
            state,_=spice.spkez(target,et,'J2000','NONE',center)
            rows.append([target,center,format(et,'.17g'),*[format(v,'.17g') for v in state]])
    dest=Path(__file__).resolve().parent/'reference'/'spk_de440s.csv'
    tmp=None
    try:
        with tempfile.NamedTemporaryFile('w',dir=dest.parent,delete=False,newline='') as out:
            tmp=Path(out.name)
            out.write('# CSPICE '+spice.tkvrsn('TOOLKIT')+' via SpiceyPy '+spice.__version__+'; NONE, J2000/ICRF.\n')
            out.write('# DE440s SHA256 '+digest+'\n')
            out.write('# target,center,TDB_seconds_from_J2000,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s\n')
            csv.writer(out,lineterminator='\n').writerows(rows)
        tmp.replace(dest)
    finally:
        if tmp is not None:tmp.unlink(missing_ok=True)
        spice.kclear()
    print(f'{dest}: {len(rows)} independent states')


if __name__=='__main__':main()
