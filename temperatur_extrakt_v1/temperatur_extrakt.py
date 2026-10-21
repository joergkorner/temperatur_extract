#!/usr/bin/env python3
"""
temperatur_extrakt.py — pulls out of every IGC flight log what a hypsometric temperature analysis needs, and nothing else.

Companion of flightstates.py (the segmenter), run on the same archive: the segmenter keeps one altitude and drops the
second, this script keeps what the segmenter drops. It writes NO track: no positions beyond a 1-km mean per layer,
no file names unless asked, no pilot. Two tables come out:

  PREFIX_fluege.csv.gz     one row per log (also logs with a single altitude, so the archive is described completely)
      id, date, t_first, lat_first, lon_first, n_fixes, dt_median_s, two_altitudes, recorder, app_version, family,
      d_median_m, d_mad_m, z_min, z_max
        id            your key (--id hash = 10 hex digits of the file name, --id name = file name, --id TEXT = constant)
        t_first       second of the day (UTC) of the first B record; lat_first, lon_first its position x 1e5
                      (join key to the segmenter's state line, whose first segment starts at or shortly after it)
        two_altitudes 1 if the log carries a pressure altitude distinct from the GNSS altitude (same rule as the paper)
        recorder      HFFTY record of the header, 30 characters; app_version: HFRFW record, 20 characters
        family        instrument family from the recorder string (Flytec, Flymaster, Syride, XC Tracer, Skytraxx,
                      Naviter, AIR3, Compass, Skydrop, Skybean, LXNav, Digifly, Phone, other, ?)
        d_median_m    median of (pressure altitude - GNSS altitude) over the flight; d_mad_m its scaled MAD
                      (a log whose pressure altitude is tied to GNSS shows |d_median| < 30 m)
  PREFIX_schichten.csv.gz  one row per log, time window and 25-hPa layer (logs with two altitudes only)
      id, date, window, p_hPa, Tv_C, n_fixes, p_span_hPa, lat, lon, z_m, extent_km, hour_utc, share_climbing, recalibrating
        window        7-10, 10-14, 14-18 (broad windows, UTC) and h06 ... h19 (single hours)
        p_hPa         layer centre; layers of 25 hPa between 875 and 525 hPa
        Tv_C          mean virtual temperature of the layer from the regression of GNSS altitude on ln p (hypsometric
                      equation), from at least 40 fixes spanning at least 15 hPa
        lat, lon      mean position of the fixes in the layer, rounded to 0.01 deg (about 1 km)
        share_climbing share of fixes with barometric vertical speed (15-s mean) above +0.3 m/s
        recalibrating 1 if the layers of this flight-window reproduce the standard atmosphere (|median difference| < 1.5 K
                      over at least two layers) — an instrument that adjusts its barometer to GNSS; drop these rows

Method and thresholds: Korner (2026), Paraglider flight logs as hypsometers, https://doi.org/10.31223/X5PJ7S.
Dependencies: numpy only. Python 3.8+.

    python3 temperatur_extrakt.py IGC_FOLDER PREFIX [--kerne 8] [--id hash|name|TEXT]
    python3 temperatur_extrakt.py IGC_FOLDER PREFIX --probe 500      (first 500 files only, for a test)

Files are found recursively (*.igc, *.IGC). A corrupt file is skipped with a note on stderr. Progress every 10 000 files.
Author: Jörg Korner, 2026. Free to use.
"""
import sys, os, re, gzip, time, hashlib, glob
import numpy as np
from multiprocessing import Pool

R_D, G = 287.05, 9.80665
LAYER_EDGES = np.arange(875, 524, -25)
WINDOWS = {"7-10": (7, 10), "10-14": (10, 14), "14-18": (14, 18)}
WINDOWS.update({f"h{h:02d}": (h, h + 1) for h in range(6, 20)})
MIN_FIXES_FLIGHT, MIN_FIXES_WINDOW, MIN_FIXES_LAYER, MIN_SPAN_HPA = 300, 100, 40, 15
W_CLIMB = 0.3; RECAL_K = 1.5; MAX_FILE_BYTES = 30_000_000; MAX_ROLLOVERS = 3
PHONES = ("samsung", "pixel", "google", "client:", "flyskyhy", "iphone", "apple", "xiaomi", "oneplus", "tcl", "huawei",
          "motorola", "nokia", "hmd", "oppo", "sony", "fairphone", "realme", "vivo", "honor")
FAMILIES = (("xc_tracer", "XC Tracer"), ("xctracer", "XC Tracer"), ("air3", "AIR3"), ("flytec", "Flytec"), ("flymaster", "Flymaster"),
            ("syride", "Syride"), ("naviter", "Naviter"), ("oudie", "Naviter"), ("skytraxx", "Skytraxx"), ("brauniger", "Flytec"),
            ("bräuniger", "Flytec"), ("compass", "Compass"), ("skydrop", "Skydrop"), ("skybean", "Skybean"), ("lxnav", "LXNav"), ("digifly", "Digifly"))

def family(recorder, fw):
    s = (recorder + " " + fw).lower()
    if not recorder.strip() or recorder.strip() == "?": return "?"
    for k, name in FAMILIES:
        if k in s: return name
    if any(k in s for k in PHONES) or fw.startswith("0."): return "Phone"
    return "other"

def pressure_hpa(pa): return 1013.25 * (1 - 2.25577e-5 * pa) ** 5.25588
def standard_temperature(p): return 15 - 0.0065 * 44330.8 * (1 - (p / 1013.25) ** 0.190263)
def layer_temperature(lnp, z):
    slope = np.linalg.lstsq(np.vstack([lnp, np.ones(len(lnp))]).T, z, rcond=None)[0][0]
    return -slope * G / R_D - 273.15

def read_igc(path):
    if os.path.getsize(path) > MAX_FILE_BYTES: raise ValueError("file too large")
    T, LA, LO, PA, GA = [], [], [], [], []; date = ""; rec = "?"; fw = ""; nb = ndiff = nzero = 0
    with open(path, "r", encoding="latin-1", errors="ignore") as fh:
        for z in fh:
            if z.startswith("B") and len(z) >= 35:
                try:
                    t = int(z[1:3]) * 3600 + int(z[3:5]) * 60 + int(z[5:7])
                    la = (int(z[7:9]) + int(z[9:14]) / 60000.0) * (1 if z[14] == "N" else -1)
                    lo = (int(z[15:18]) + int(z[18:23]) / 60000.0) * (1 if z[23] == "E" else -1)
                    pa = int(z[25:30]); ga = int(z[30:35])
                except ValueError:
                    continue
                if nb < 400:
                    nb += 1
                    if pa == 0 or ga == 0: nzero += 1
                    elif pa != ga: ndiff += 1
                T.append(t); LA.append(la); LO.append(lo); PA.append(pa); GA.append(ga)
            elif z.startswith("HFDTE") and not date:
                d = re.sub(r"[^0-9]", "", z[5:])
                if len(d) >= 6: y = int(d[4:6]); date = f"{2000 + y if y < 80 else 1900 + y}-{d[2:4]}-{d[0:2]}"
            elif z.startswith("HFFTY"):
                rec = z.strip().split(":", 1)[-1].strip()[:30].replace(",", " ").replace(";", " ")
            elif z.startswith("HFRFW"):
                fw = z.strip().split(":", 1)[-1].strip()[:20].replace(",", " ").replace(";", " ")
    two = 1 if (nb > 0 and nzero < 0.9 * nb and ndiff >= 0.5 * nb) else 0
    return date, rec, fw, two, [np.array(x, float) for x in (T, LA, LO, PA, GA)]

def one(args):
    path, idmode = args
    name = os.path.basename(path)
    fid = hashlib.sha1(name.encode()).hexdigest()[:10] if idmode == "hash" else (name if idmode == "name" else idmode)
    try:
        date, rec, fw, two, (T, LA, LO, PA, GA) = read_igc(path)
        if len(T) < 5: return None
        back = np.diff(T) < 0
        if back.sum() > MAX_ROLLOVERS: raise ValueError("timestamps run backwards")
        T = T + np.concatenate(([0], np.cumsum(back) * 86400.0))
        dt = float(np.median(np.diff(T))) if len(T) > 1 else 0.0
        flug = [fid, date, int(T[0]), int(round(LA[0] * 1e5)), int(round(LO[0] * 1e5)), len(T), round(dt, 1), two, rec, fw, family(rec, fw)]
        rows = []
        if not two:
            return flug + ["", "", "", ""], rows
        ok = (PA != 0) & (GA != 0); T, LA, LO, PA, GA = T[ok], LA[ok], LO[ok], PA[ok], GA[ok]
        keep = np.concatenate([[True], np.diff(T) > 0]); T, LA, LO, PA, GA = T[keep], LA[keep], LO[keep], PA[keep], GA[keep]
        d = PA - GA; dm = float(np.median(d)); dmad = float(1.4826 * np.median(np.abs(d - dm)))
        flug += [round(dm, 1), round(dmad, 1), int(PA.min()), int(PA.max())]
        if len(T) < MIN_FIXES_FLIGHT: return flug, rows
        p = pressure_hpa(PA); lnp = np.log(p)
        w = np.convolve(np.gradient(PA, T), np.ones(15) / 15, mode="same")
        x = np.cos(np.radians(LA)) * np.radians(LO) * 6371000; y = np.radians(LA) * 6371000
        for win, (h0, h1) in WINDOWS.items():
            m = (T >= h0 * 3600) & (T < h1 * 3600)
            if m.sum() < MIN_FIXES_WINDOW: continue
            vals = []
            for k in range(len(LAYER_EDGES) - 1):
                hi, lo = LAYER_EDGES[k], LAYER_EDGES[k + 1]; s = m & (p <= hi) & (p > lo)
                n = int(s.sum())
                if n < MIN_FIXES_LAYER or (p[s].max() - p[s].min()) < MIN_SPAN_HPA: continue
                Tv = layer_temperature(lnp[s], GA[s]); pc = (hi + lo) / 2
                ext = float(np.hypot(x[s].max() - x[s].min(), y[s].max() - y[s].min())) / 1000
                vals.append([fid, date, win, pc, round(Tv, 2), n, round(float(p[s].max() - p[s].min()), 1), round(float(LA[s].mean()), 2),
                             round(float(LO[s].mean()), 2), int(GA[s].mean()), round(ext, 1), round(float(T[s].mean() / 3600), 2),
                             round(float(np.mean(w[s] > W_CLIMB)), 2), Tv - standard_temperature(pc)])
            if not vals: continue
            recal = 1 if (len(vals) >= 2 and abs(float(np.median([v[-1] for v in vals]))) < RECAL_K) else 0
            for v in vals: rows.append(v[:-1] + [recal])
        return flug, rows
    except Exception as e:
        sys.stderr.write(f"skipped {name}: {e}\n"); return None

def main():
    a = sys.argv[1:]
    if len(a) < 2: print(__doc__); sys.exit(1)
    folder, prefix = a[0], a[1]
    kerne = int(a[a.index("--kerne") + 1]) if "--kerne" in a else max(1, os.cpu_count() - 1)
    idmode = a[a.index("--id") + 1] if "--id" in a else "hash"
    probe = int(a[a.index("--probe") + 1]) if "--probe" in a else 0
    files = sorted(glob.glob(os.path.join(folder, "**", "*.igc"), recursive=True) + glob.glob(os.path.join(folder, "**", "*.IGC"), recursive=True))
    if probe: files = files[:probe]
    print(f"{len(files)} files, {kerne} cores", flush=True); t0 = time.time()
    fh_f = gzip.open(prefix + "_fluege.csv.gz", "wt"); fh_s = gzip.open(prefix + "_schichten.csv.gz", "wt")
    fh_f.write("id,date,t_first,lat_first,lon_first,n_fixes,dt_median_s,two_altitudes,recorder,app_version,family,d_median_m,d_mad_m,z_min,z_max\n")
    fh_s.write("id,date,window,p_hPa,Tv_C,n_fixes,p_span_hPa,lat,lon,z_m,extent_km,hour_utc,share_climbing,recalibrating\n")
    nf = ns = 0
    with Pool(kerne) as pool:
        for i, r in enumerate(pool.imap_unordered(one, [(f, idmode) for f in files], chunksize=50)):
            if r is None: continue
            flug, rows = r; nf += 1; ns += len(rows)
            fh_f.write(",".join(str(v) for v in flug) + "\n")
            for v in rows: fh_s.write(",".join(str(x) for x in v) + "\n")
            if i % 10000 == 0: print(f"{i} files, {ns} layer rows, {round(time.time() - t0)} s", flush=True)
    fh_f.close(); fh_s.close()
    print(f"done: {nf} flights, {ns} layer rows, {round(time.time() - t0)} s -> {prefix}_fluege.csv.gz, {prefix}_schichten.csv.gz")

if __name__ == "__main__":
    main()
