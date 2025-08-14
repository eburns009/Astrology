# -*- coding: utf-8 -*-
"""
Astrology Chart App — Fagan/Bradley (SVP) Sidereal
- Flask + Skyfield
- True geometric, geocentric longitudes
- Fagan/Bradley ayanamsa (precise fallback): J1950 base + rate
- HTML lives in templates/chart.html
"""

from __future__ import annotations
from datetime import datetime, timezone as _tz, timedelta
from functools import lru_cache
from math import floor
from pathlib import Path

import pytz
from flask import Flask, render_template, request
from skyfield.api import Loader
from skyfield.framelib import ecliptic_frame as ECLIPTIC_FRAME

app = Flask(__name__)

# ---------------- Skyfield cache ----------------
DATA_DIR = Path.home() / ".skyfield-data"
_loader = Loader(str(DATA_DIR))

@lru_cache(maxsize=1)
def get_timescale():
    return _loader.timescale()

@lru_cache(maxsize=1)
def get_ephemeris():
    # small & accurate; change to de440s if you prefer
    return _loader("de421.bsp")

# ---------------- Utilities ----------------
ZODIAC_SIGNS = [
    "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
    "Libra","Scorpio","Sagittarius","Capricorn","Aquarius","Pisces"
]

PLANETS = [
    ("Sun", "sun"),
    ("Moon", "moon"),
    ("Mercury", "mercury"),
    ("Venus", "venus"),
    ("Earth", "earth"),
    ("Mars", "mars barycenter"),
    ("Jupiter", "jupiter barycenter"),
    ("Saturn", "saturn barycenter"),
    ("Uranus", "uranus barycenter"),
    ("Neptune", "neptune barycenter"),
    ("Pluto", "pluto barycenter"),
]

def normalize_deg(x: float) -> float:
    x = x % 360.0
    return x if x >= 0.0 else x + 360.0

def format_longitude(lon: float) -> str:
    lon = normalize_deg(lon)
    sign_idx = int(lon // 30)
    sign = ZODIAC_SIGNS[sign_idx]
    in_sign = lon - sign_idx * 30
    deg = int(floor(in_sign))
    minutes_full = (in_sign - deg) * 60.0
    minutes = int(floor(minutes_full))
    seconds = int(round((minutes_full - minutes) * 60.0))
    return f"{deg:02d}° {minutes:02d}' {seconds:02d}\" {sign}"

def julian_day_utc(dt: datetime) -> float:
    """Meeus JD for a UTC datetime."""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=_tz.utc)
    else:
        dt = dt.astimezone(_tz.utc)
    y, m = dt.year, dt.month
    D = dt.day + (dt.hour + (dt.minute + dt.second/60.0)/60.0) / 24.0
    if m <= 2:
        y -= 1; m += 12
    if (y > 1582) or (y == 1582 and m > 10) or (y == 1582 and m == 10 and D >= 15):
        A = y // 100
        B = 2 - A + A // 4
    else:
        B = 0
    return int(365.25*(y+4716)) + int(30.6001*(m+1)) + D + B - 1524.5

def fagan_bradley_ayanamsa(dt: datetime) -> float:
    """
    Fagan/Bradley (SVP) precise fallback:
      - Base at J1950.0 = 24°02′39.43″ = 24.044286111°
      - Rate = 50.290966″/yr
      - Year length = 365.24219879 days (tropical)
    """
    dt_utc = dt if dt.tzinfo and dt.utcoffset() == timedelta(0) else dt.astimezone(_tz.utc)
    jd = julian_day_utc(dt_utc)
    J1950 = 2433282.5
    years = (jd - J1950) / 365.24219879
    base = 24.044286111
    rate = 50.290966 / 3600.0  # arcsec/yr -> deg/yr
    return normalize_deg(base + years * rate)

def planetary_longitudes(dt_utc: datetime, helio: bool = False) -> dict:
    """True geometric ecliptic longitudes, geocentric by default (no topocentric)."""
    ts = get_timescale()
    t = ts.from_datetime(dt_utc if dt_utc.tzinfo else dt_utc.replace(tzinfo=_tz.utc))
    eph = get_ephemeris()
    origin = eph["sun"] if helio else eph["earth"]
    out = {}
    for name, key in PLANETS:
        if key == "earth" and not helio:
            continue
        target = eph[key]
        astrometric = origin.at(t).observe(target)  # geometric (no aberration/deflection)
        lat_ecl, lon_ecl, _ = astrometric.frame_latlon(ECLIPTIC_FRAME)
        out[name] = normalize_deg(lon_ecl.degrees)
    return out

# ---------------- Route ----------------
@app.route("/", methods=["GET", "POST"])
def chart():
    error = None
    data = None

    # defaults shown in the form
    defaults = {
        "person": "",
        "place": "Fort Knox, KY",
        "dt": "1962-07-02T23:33",
        "tz": "America/New_York",
        "lat": "37.897655",
        "lon": "-85.903852",
        "zodiac": "sidereal",   # default to Fagan/Bradley
        "helio": "0",
    }

    if request.method == "POST":
        try:
            person = request.form.get("person", "").strip()
            place = request.form.get("place","").strip()
            dt_str = request.form.get("dt", "").strip()               # YYYY-MM-DDTHH:MM
            tz_name = request.form.get("tz", defaults["tz"]).strip()
            lat = float(request.form.get("lat", defaults["lat"]))
            lon = float(request.form.get("lon", defaults["lon"]))
            zodiac = request.form.get("zodiac", "sidereal")           # "sidereal" or "tropical"
            helio = request.form.get("helio", "0")                    # "0" or "1"

            if not dt_str:
                raise ValueError("Please provide a local birth date & time.")

            tz = pytz.timezone(tz_name)
            local_dt = datetime.strptime(dt_str, "%Y-%m-%dT%H:%M")
            local_dt = tz.localize(local_dt, is_dst=None)
            dt_utc = local_dt.astimezone(pytz.UTC).replace(tzinfo=_tz.utc)

            longs_trop = planetary_longitudes(dt_utc, helio=(helio == "1"))
            if zodiac == "sidereal":
                ay = fagan_bradley_ayanamsa(dt_utc)  # Fagan/Bradley only
                longs = {k: normalize_deg(v - ay) for k, v in longs_trop.items()}
                ay_fmt = f"{ay:.6f}°"
                zodiac_label = "Sidereal (Fagan/Bradley)"
            else:
                longs = longs_trop
                ay_fmt = "—"
                zodiac_label = "Tropical"

            rows = []
            for name, _key in PLANETS:
                if name == "Earth" and helio != "1":
                    continue
                if name not in longs:
                    continue
                lon_deg = longs[name]
                rows.append({
                    "name": name,
                    "lon_fmt": format_longitude(lon_deg),
                    "sign": ZODIAC_SIGNS[int(lon_deg // 30)],
                    "decimal_degrees": round(lon_deg, 5),
                })

            data = {
                "person": person or "—",
                "place": place or "—",
                "dt_disp": local_dt.strftime("%Y-%m-%d %H:%M"),
                "tz": tz_name,
                "utc_offset": local_dt.utcoffset().total_seconds() / 3600.0,
                "dst": "Yes" if local_dt.dst() and local_dt.dst().total_seconds() != 0 else "No",
                "lst": "—",
                "lat": round(lat, 6),
                "lon": round(lon, 6),
                "frame": "Heliocentric" if helio == "1" else "Geocentric",
                "zodiac": zodiac_label,
                "ayanamsa": ay_fmt,
                "table": rows,
            }

        except Exception as e:
            error = str(e)

    return render_template("chart.html",
                           data=data,
                           error=error,
                           defaults=defaults)

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 10000))
    app.run(host="0.0.0.0", port=port)
