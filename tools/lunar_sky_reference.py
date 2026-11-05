#!/usr/bin/env python3
"""Regenerate independent JPL/USNO lunar fixtures; never evaluate Spp.

python3 tools/lunar_sky_reference.py [--cache /path/to/downloaded/responses]
Normal builds/tests use committed CSVs and need no network or Python packages.
"""
import argparse
import csv
from contextlib import contextmanager
import datetime
import json
from pathlib import Path
import time
import tempfile
import urllib.parse
import urllib.request

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cache', type=Path, required=True,
                    help='cache directory for downloaded, unmodified API responses')
args = parser.parse_args()
args.cache.mkdir(parents=True, exist_ok=True)
dest = Path(__file__).resolve().parent / 'reference'
dates = ['2000-01-01 12:00', '2026-11-03 12:00', '2026-11-04 00:00',
         '2024-04-08 18:00', '2026-11-24 12:00', '2050-06-01 00:00',
         '1900-01-01 12:00', '2100-06-01 00:00']
jd_list = ' '.join(str(datetime.datetime.fromisoformat(d).replace(
    tzinfo=datetime.timezone.utc).timestamp()/86400+2440587.5) for d in dates)

def fetch(name, url):
    path = args.cache / name
    if not path.exists():
        for attempt in range(4):
            try:
                content = urllib.request.urlopen(url, timeout=60).read()
                value = json.loads(content)
                if value.get('error'):
                    raise ValueError(value['error'])
                path.write_bytes(content)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(1)
    value = json.loads(path.read_bytes())
    if value.get('error'):
        raise ValueError(f'{path}: {value["error"]}')
    return value

@contextmanager
def fixture(name):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile('w', newline='', dir=dest, suffix='.tmp', delete=False) as file:
            temporary = Path(file.name)
            yield file
        temporary.replace(dest / name)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

with fixture('moon_horizons.csv') as file:
    file.write('# JPL Horizons observer tables, DE441; AIRLESS, UTC, apparent true-of-date.\n'
               '# Modern epochs only: TT-UTC and DUT1 conventions match Spp.\n'
               '# source,jd_utc,lat,lon,elev_m,dut1,ra,dec,az,alt,range_km,diameter_arcsec,illum_fraction,elongation,phase_angle_phi,bright_limb_deg\n')
    writer = csv.writer(file, lineterminator='\n')
    for name, lon, lat, height in [('paris', 2.3522, 48.8566, 0.035),
                                  ('south', 151.2093, -33.8688, 0.058),
                                  ('polar', -156.7886, 71.2906, 0.0),
                                  ('geocenter', 0, 0, -6378.137)]:
        params = dict(format='json', COMMAND="'301'", EPHEM_TYPE="'OBSERVER'",
            CENTER="'500@399'" if name == 'geocenter' else "'coord@399'",
            COORD_TYPE="'GEODETIC'", SITE_COORD=f"'{lon},{lat},{height}'",
            TLIST=f"'{jd_list}'", TIME_TYPE="'UT'", CAL_FORMAT="'BOTH'",
            ANG_FORMAT="'DEG'", APPARENT="'AIRLESS'", EXTRA_PREC="'YES'",
            REF_SYSTEM="'ICRF'",
            QUANTITIES="'2,4,7,10,13,16,20,23,24,31,43,49'", CSV_FORMAT="'YES'", OBJ_DATA="'NO'")
        data = fetch(name+'-disk43.json', 'https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(params))
        result = data['result']
        if 'DE441' not in result:
            raise ValueError('Unexpected Horizons ephemeris; review before regenerating')
        rows = csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines())
        for row in rows:
            row = [v.strip() for v in row]
            jd = float(row[1])
            if not 2451544.5 <= jd < 2461406.5:
                continue
            val = lambda col: float(row[col]) if row[col] != 'n.a.' else 'nan'
            writer.writerow([name, jd, lat, lon, height*1000, val(23), val(4), val(5),
                             val(6), val(7), val(13)*149597870.7, val(10), val(9)/100,
                             val(15), val(20), val(11)])

with fixture('moon_usno_phases.csv') as file:
    file.write('# USNO api/moon/phases/year?year=YEAR. Times are Universal Time, at minute precision.\n'
               '# year,month,day,hour,minute,phase (0=new,1=first quarter,2=full,3=last quarter)\n')
    writer = csv.writer(file, lineterminator='\n')
    kinds = {'New Moon': 0, 'First Quarter': 1, 'Full Moon': 2, 'Last Quarter': 3}
    for year in [2000, 2024, 2026, 2050]:
        data = fetch(f'phases-{year}.json', f'https://aa.usno.navy.mil/api/moon/phases/year?year={year}')
        for event in data['phasedata']:
            h, m = map(int, event['time'].split(':'))
            writer.writerow([event['year'], event['month'], event['day'], h, m, kinds[event['phase']]])

print('Wrote moon_horizons.csv and moon_usno_phases.csv')
