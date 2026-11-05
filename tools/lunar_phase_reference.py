#!/usr/bin/env python3
"""Generate matched-TT JPL phase short arcs, independently of Spp.

This explicitly requested regeneration downloads two observer tables. Normal
builds and benchmarks use the committed CSV without Python or network access.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    reference = Path(__file__).resolve().parent / 'reference'
    cases, counts = [], {}
    with (reference / 'moon_usno_phases.csv').open() as source:
        for row in csv.reader(line for line in source if not line.startswith('#')):
            year, month, day, hour, minute, kind = map(int, row)
            counts[year] = counts.get(year, 0) + 1
            if year != 2026 and counts[year] > 4:
                continue
            # USNO only supplies the neighborhood. TT, not UT/UTC, is the
            # independent table's time scale; no leap/Delta-T model is shared.
            center = (dt.datetime(year, month, day, hour, minute, tzinfo=dt.timezone.utc)
                      .timestamp() / 86400 + 2440587.5 + 90 / 86400)
            cases.append((kind, [center + seconds / 86400 for seconds in (-180, 0, 90, 180)]))
    epochs = [jd for _, arc in cases for jd in arc]

    def batch(body, samples):
        params = dict(format='json', COMMAND=f"'{body}'", OBJ_DATA="'NO'",
                      EPHEM_TYPE="'OBSERVER'", CENTER="'500@399'",
                      TIME_TYPE="'TT'", TLIST_TYPE="'JD'",
                      TLIST="'" + ' '.join(f'{jd:.12f}' for jd in samples) + "'",
                      QUANTITIES="'31'", APPARENT="'AIRLESS'", REF_SYSTEM="'ICRF'",
                      ANG_FORMAT="'DEG'", EXTRA_PREC="'YES'", CSV_FORMAT="'YES'",
                      CAL_FORMAT="'BOTH'", TIME_DIGITS="'FRACSEC'")
        url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(params)
        digest = hashlib.sha256(url.encode()).hexdigest()[:16]
        path = args.cache / f'phase-tt-{body}-{digest}.json'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                content = response.read()
            value = json.loads(content)
            if value.get('error'):
                raise ValueError(value['error'])
            path.write_bytes(content)
        value = json.loads(path.read_bytes())
        if value.get('error'):
            raise ValueError(value['error'])
        result = value['result']
        if 'DE441' not in result or 'Date__(TT)' not in result:
            raise ValueError(f'{path}: review ephemeris/time-scale change before accepting')
        rows = list(csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))
        if len(rows) != len(samples):
            raise ValueError('Horizons epoch count differs from the request')
        return [(float(row[1]), float(row[4])) for row in rows]

    def table(body):
        rows = []
        for offset in range(0, len(epochs), 32):
            rows.extend(batch(body, epochs[offset:offset+32]))
        return rows

    moon, sun = table(301), table(10)
    path = reference / 'moon_horizons_phase_tt.csv'
    temporary = None
    try:
        with tempfile.NamedTemporaryFile('w', newline='', dir=reference, suffix='.tmp', delete=False) as output:
            temporary = Path(output.name)
            output.write('# JPL Horizons DE441, geocenter, airless apparent IAU76/80 ecliptic of date.\n'
                         '# TT input/output; four samples per event; USNO only locates the neighborhood.\n'
                         '# case,phase,jd_tt,moon_longitude_deg,sun_longitude_deg\n')
            writer = csv.writer(output, lineterminator='\n')
            for i, (kind, arc) in enumerate(cases):
                for j in range(4):
                    m, s = moon[4*i+j], sun[4*i+j]
                    if abs(m[0]-arc[j]) > 1e-8 or m[0] != s[0]:
                        raise ValueError('Horizons epoch mismatch')
                    writer.writerow((i, kind, f'{m[0]:.12f}', f'{m[1]:.12f}', f'{s[1]:.12f}'))
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(f'Wrote {path.name}: {len(cases)} independent phase neighborhoods')


if __name__ == '__main__':
    main()
