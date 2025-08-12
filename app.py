"""
New Astrology Emerging — Switchboard (Local) with Houses + Aspects + Birthplace
Now supports:
- Birth time (user-entered local time)
- Birthplace by city (auto geocoding to lat/lon)
- Auto timezone detection from birthplace (or manual choose)
- Text listing for planets and aspects (in addition to tables)
"""

from __future__ import annotations
from math import atan2, degrees, radians, sin, cos, floor, asin
from datetime import datetime, timezone as _tz
import pytz

from flask import Flask, request, render_template_string, jsonify

# Skyfield
from skyfield.api import load, wgs84, Loader
from skyfield.framelib import ecliptic_frame as ECLIPTIC_FRAME  # fixed ecliptic frame access

# Geocoding + timezone lookup
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder

# Initialize helpers
app = Flask(__name__)
_geocoder = Nominatim(user_agent="nae-switchboard")
_tzf = TimezoneFinder()

_loader = None
_ts = None
_eph = None

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

ZODIAC_SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"
]

ASPECTS_DEF = [
    {"name":"Conjunction","angle":0,   "key":"conj", "default_orb":8},
    {"name":"Opposition", "angle":180, "key":"opp",  "default_orb":8},
    {"name":"Trine",      "angle":120, "key":"tri",  "default_orb":6},
    {"name":"Square",     "angle":90,  "key":"sqr",  "default_orb":6},
    {"name":"Sextile",    "angle":60,  "key":"sex",  "default_orb":4},
    {"name":"Quincunx",   "angle":150, "key":"qun",  "default_orb":3},
]

# -------- utilities --------
def normalize_deg(angle: float) -> float:
    a = angle % 360.0
    return a if a >= 0 else a + 360.0

def dms(angle: float) -> tuple[int, int, float]:
    a = normalize_deg(angle)
    d = floor(a)
    m_full = (a - d) * 60.0
    m = floor(m_full)
    s = (m_full - m) * 60.0
    return d, m, s

def format_longitude(angle: float) -> str:
    d, m, s = dms(angle)
    sign_index = int(floor(angle / 30.0)) % 12
    in_sign = angle % 30.0
    di, mi, si = dms(in_sign)
    return f"{ZODIAC_SIGNS[sign_index]} {di:02d}°{mi:02d}'{si:04.1f}″"

def _is_leap(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

def fagan_bradley_ayanamsa(dt: datetime) -> float:
    """High-accuracy Fagan/Bradley ayanamsha.
    Tries Swiss Ephemeris (pyswisseph) first; falls back to a smooth
    polynomial approximation if the module isn't available.
    Returns degrees in [0,360).
    """
    # Work in UTC
    if dt.tzinfo is None:
        dt_utc = dt.replace(tzinfo=_tz.utc)
    else:
        dt_utc = dt.astimezone(_tz.utc)

    # Preferred: Swiss Ephemeris (very close to canonical values)
    try:
        import swisseph as swe  # provided by the 'pyswisseph' package
        jd = swe.julday(
            dt_utc.year,
            dt_utc.month,
            dt_utc.day,
            dt_utc.hour + dt_utc.minute/60.0 + dt_utc.second/3600.0,
            swe.GREG_CAL,
        )
        swe.set_sid_mode(swe.SIDM_FAGAN_BRADLEY, 0, 0)
        ay = float(swe.get_ayanamsa_ut(jd))
        return normalize_deg(ay)
    except Exception:
        # Fallback: smooth polynomial around J2000 (good but not SE-precise)
        y = dt_utc.year + (
            dt_utc.timetuple().tm_yday - 1 +
            (dt_utc.hour + dt_utc.minute/60 + dt_utc.second/3600)/24
        ) / (366 if _is_leap(dt_utc.year) else 365)
        t = y - 2000.0
        base = 24.042044444
        rate = -0.013004167
        curve = -0.000000164 * t * t
        ay = base + rate * t + curve
        return normalize_deg(ay)

# -------- Skyfield helpers --------
def get_loader() -> Loader:
    global _loader
    if _loader is None:
        _loader = load
    return _loader

def get_timescale():
    global _ts
    if _ts is None:
        _ts = get_loader().timescale()
    return _ts

def get_ephemeris():
    global _eph
    if _eph is None:
        _eph = get_loader()("de440s.bsp")
    return _eph

# -------- core calculations --------
def planetary_longitudes(dt: datetime, lat: float, lon: float, elevation_m: float,
                         helio: bool = False) -> dict:
    ts = get_timescale()
    t = ts.from_datetime(dt)
    eph = get_ephemeris()

    if helio:
        origin = eph["sun"]
    else:
        topos = wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=elevation_m)
        origin = eph["earth"] + topos   # Earth + observer (fix for .observe())

    results = {}
    for name, key in PLANETS:
        if key == "earth":
            if not helio:
                continue
        target = eph[key]
        astrometric = origin.at(t).observe(target)
        lat_ecl, lon_ecl, distance = astrometric.frame_latlon(ECLIPTIC_FRAME)
        lam = normalize_deg(lon_ecl.degrees)
        results[name] = lam
    return results

def to_sidereal(longitudes_tropical: dict, dt: datetime) -> dict:
    ay = fagan_bradley_ayanamsa(dt)
    return {name: normalize_deg(lon - ay) for name, lon in longitudes_tropical.items()}

# -------- ASC/MC/Houses --------
# High-precision obliquity utilities (fallback to Laskar mean obliquity)
OBLIQ_DEG_MEAN_J2000 = 23.43929111  # fallback constant

def _julian_day_utc(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=_tz.utc)
    else:
        dt = dt.astimezone(_tz.utc)
    y, m = dt.year, dt.month
    D = dt.day + (dt.hour + dt.minute/60 + dt.second/3600)/24.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 5
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + D + B - 1524.5
    return JD

def mean_obliquity_laskar(dt: datetime) -> float:
    T = (_julian_day_utc(dt) - 2451545.0) / 36525.0
    seconds = (84381.406
               - 46.836769*T
               - 0.0001831*(T**2)
               + 0.00200340*(T**3)
               - 5.76e-7*(T**4)
               - 4.34e-8*(T**5))
    return seconds / 3600.0  # degrees

def obliquity_deg(dt: datetime) -> float:
    """Return obliquity of the ecliptic in degrees.
    Tries Swiss Ephemeris if available; otherwise uses Laskar mean obliquity.
    """
    try:
        import swisseph as swe  # type: ignore
        jd = _julian_day_utc(dt)
        # Swiss provides mean obliquity with swe.obl_ecl(jd, flag)
        # flag = 0 -> mean obliquity of date; flag = 1 -> true obliquity (includes nutation)
        eps, _, _ = swe.obl_ecl(jd, 1)
        return float(eps)
    except Exception:
        return mean_obliquity_laskar(dt)

def gmst_hours(dt: datetime) -> float:
    t = get_timescale().from_datetime(dt)
    return t.gast  # hours (apparent sidereal time incl. nutation) 

def lst_deg(dt: datetime, lon_deg: float) -> float:
    lst_h = gmst_hours(dt) + (lon_deg / 15.0)
    lst_h %= 24.0
    return lst_h * 15.0

def asc_mc(dt: datetime, lat_deg: float, lon_deg: float) -> tuple[float, float]:
    """Return (ASC, MC) ecliptic longitudes in degrees.
    Tries Swiss Ephemeris for top-precision. Falls back to analytic
    method using apparent sidereal time and obliquity-of-date.
    """
    # 1) Try Swiss Ephemeris (houses_ex with Equal houses just to get ASC/MC)
    try:
        import swisseph as swe  # type: ignore
        jd = _julian_day_utc(dt)
        flags = 0
        cusps, ascmc = swe.houses_ex(jd, flags, float(lat_deg), float(lon_deg), b'E')
        asc = float(ascmc[0])  # Ascendant
        mc = float(ascmc[1])   # Midheaven
        return normalize_deg(asc), normalize_deg(mc)
    except Exception:
        pass

    # 2) Fallback analytic method
    ε = radians(obliquity_deg(dt))
    φ = radians(lat_deg)
    θ = radians(lst_deg(dt, lon_deg))

    # MC (analytic)
    lam_mc = degrees(atan2(sin(θ), cos(θ)*cos(ε)))
    lam_mc = normalize_deg(lam_mc)

    # ASC via sampling/root-refinement (altitude ≈ 0°, eastern)
    def eq_to_altaz(alpha, delta):
        H = (θ - alpha)
        while H > 3.141592653589793: H -= 2*3.141592653589793
        while H < -3.141592653589793: H += 2*3.141592653589793
        sin_h = sin(φ)*sin(delta) + cos(φ)*cos(delta)*cos(H)
        h = asin(max(-1.0, min(1.0, sin_h)))
        sinA = -sin(H) * cos(delta)
        cosA = (sin(delta)) / max(1e-9, cos(φ))
        A = atan2(sinA, cosA)
        if A < 0: A += 2*3.141592653589793
        return h, A

    def ecl_to_eq(lam):
        lamr = radians(lam)
        alpha = atan2(sin(lamr)*cos(ε), cos(lamr))
        if alpha < 0: alpha += 2*3.141592653589793
        delta = asin(sin(lamr)*sin(ε))
        return alpha, delta

    prev_h = None
    prev_lam = None
    roots = []
    for lam in [i*1.0 for i in range(361)]:
        alpha, delta = ecl_to_eq(lam)
        h, A = eq_to_altaz(alpha, delta)
        if prev_h is not None and (h*prev_h) < 0:
            lam0 = prev_lam + (0 - prev_h) * (lam - prev_lam) / (h - prev_h)
            alpha0, delta0 = ecl_to_eq(lam0)
            h0, A0 = eq_to_altaz(alpha0, delta0)
            roots.append((lam0, A0))
        prev_h = h
        prev_lam = lam

    asc_candidates = [normalize_deg(l) for (l, A0) in roots if 0 < degrees(A0) < 180]
    lam_asc = asc_candidates[0] if asc_candidates else 0.0
    return lam_asc, lam_mc

def equal_house_cusps(asc_deg: float, mode: str = "asc_middle") -> list[float]:
    if mode == "asc_middle":
        cusp1 = normalize_deg(asc_deg - 15.0)
    else:
        cusp1 = normalize_deg(asc_deg)
    return [normalize_deg(cusp1 + i*30.0) for i in range(12)]

# -------- Aspects --------
def angle_sep(a: float, b: float) -> float:
    d = abs(a - b) % 360.0
    return d if d <= 180 else 360 - d

def find_aspects(longs: dict, aspect_opts: dict) -> list[dict]:
    names = [n for n,_ in PLANETS if n != 'Earth']
    results = []
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            n1, n2 = names[i], names[j]
            if n1 not in longs or n2 not in longs:
                continue
            d = angle_sep(longs[n1], longs[n2])
            for spec in ASPECTS_DEF:
                if not aspect_opts.get(spec['key']+'_on', True):
                    continue
                orb = float(aspect_opts.get(spec['key']+'_orb', spec['default_orb']))
                if abs(d - spec['angle']) <= orb:
                    results.append({
                        'p1': n1, 'p2': n2,
                        'type': spec['name'], 'angle': spec['angle'], 'delta': round(d - spec['angle'], 3)
                    })
    return results

# -------- UI --------
LAYOUT = """


ABOUT = """
<!doctype html>
<html><head><meta charset="utf-8"><title>About</title>
<style>body{background:#ffffff;color:#111827;font-family:ui-sans-serif} .wrap{max-width:800px;margin:40px auto;padding:0 16px} a{color:#2563eb}</style>
</head>
<body><div class="wrap">
<h1>About & current limits</h1>
<ul>
<li><b>Accuracy:</b> Skyfield + DE440s (JPL). Swiss Ephemeris used for ayanamsha and ASC/MC if available.</li>
<li><b>Sidereal:</b> Fagan/Bradley ayanamsa uses Swiss Ephemeris when installed (fallback to polynomial).</li>
<li><b>Houses:</b> Equal only (for now). Other systems possible later.</li>
<li><b>Aspects:</b> Six standard aspects with custom orbs.</li>
<li><b>Privacy:</b> Runs locally / on your Render instance.</li>
</ul>
<p><a href="/">Back</a></p>
</div></body></html>
"""
<!doctype html>
<html><head><meta charset="utf-8"><title>About</title>
<style>body{background:#0b0f14;color:#eaf2ff;font-family:ui-sans-serif} .wrap{max-width:800px;margin:40px auto;padding:0 16px} a{color:#7cc0ff}</style>
</head>
<body><div class="wrap">
<h1>About & current limits</h1>
<ul>
<li><b>Accuracy:</b> Skyfield + DE440s (JPL). Swiss Ephemeris used for ayanamsha and ASC/MC if available.</li>
<li><b>Sidereal:</b> Fagan/Bradley ayanamsa uses a polynomial approx. Swiss Ephemeris can be integrated for exact parity.</li>
<li><b>Houses:</b> Equal only (for now). Other systems possible later.</li>
<li><b>Aspects:</b> Six standard aspects with custom orbs.</li>
<li><b>Privacy:</b> Runs locally / on your Render instance.</li>
</ul>
<p><a href="/">Back</a></p>
</div></body></html>
"""

@app.route("/about")
def about():
    return ABOUT


@app.route("/")
def index():
    """Home page with empty form and defaults."""
    default_tz = "auto"
    # Seed the datetime-local control with a valid string
    now_local = datetime.now(pytz.timezone("America/Denver")).replace(second=0, microsecond=0)
    default_dt = now_local.strftime("%Y-%m-%dT%H:%M")

    aspects = [
        {"key": spec["key"], "name": spec["name"], "orb": spec["default_orb"], "on": True}
        for spec in ASPECTS_DEF
    ]

    return render_template_string(
        LAYOUT,
        default_dt=default_dt,
        default_tz=default_tz,
        data=None,
        aspects=aspects,
    )


@app.route("/chart")
def chart():
    try:
        person = (request.args.get("person") or "").strip()
        dt_str = request.args.get("dt")
        tz_sel = request.args.get("tz", "auto")
        place = (request.args.get("place") or "").strip()
        lat_str = request.args.get("lat")
        lon_str = request.args.get("lon")
        elev_str = request.args.get("elev")
        frame = request.args.get("frame", "geo")
        zodiac = request.args.get("zodiac", "sidereal")
        house_mode = request.args.get("house_mode", "asc_middle")
        house_clockwise = request.args.get("house_clockwise", "no") == "yes"

        # Build aspect options
        aspect_opts = {}
        for spec in ASPECTS_DEF:
            k = spec["key"]
            aspect_opts[f"{k}_on"] = (request.args.get(f"{k}_on") is not None)
            try:
                aspect_opts[f"{k}_orb"] = float(request.args.get(f"{k}_orb", spec["default_orb"]))
            except Exception:
                aspect_opts[f"{k}_orb"] = spec["default_orb"]

        # Coordinates
        lat = float(lat_str) if (lat_str and lat_str.strip()) else None
        lon = float(lon_str) if (lon_str and lon_str.strip()) else None
        elev = float(elev_str) if (elev_str and elev_str.strip()) else 0.0

        if (lat is None or lon is None) and place:
            loc = _geocoder.geocode(place, addressdetails=False, language="en")
            if loc:
                lat = float(loc.latitude)
                lon = float(loc.longitude)
            else:
                # Keep them None; the template will show what we have
                pass

        if lat is None or lon is None:
            return jsonify({"error": "Please enter a valid birthplace or coordinates."}), 400

        # Timezone
        if tz_sel == "auto":
            tz_name = _tzf.timezone_at(lng=lon, lat=lat) or "UTC"
        else:
            tz_name = tz_sel
        tz = pytz.timezone(tz_name)

        # Parse local birth time and make timezone-aware
        local_dt = datetime.strptime(dt_str, "%Y-%m-%dT%H:%M")
        dt = tz.localize(local_dt)

        # Display strings
        dt_disp = dt.strftime("%Y-%m-%d %H:%M")
        offset_td = dt.utcoffset() or (dt - dt)
        total_min = int(offset_td.total_seconds() // 60)
        sign = "+" if total_min >= 0 else "-"
        hh = abs(total_min) // 60
        mm = abs(total_min) % 60
        utc_offset = f"{sign}{hh:02d}:{mm:02d}"

        # Compute positions (tropical first)
        helio = (frame == "helio")
        longs_trop = planetary_longitudes(dt, lat, lon, elev, helio=helio)

        # Sidereal conversion when requested
        ay = fagan_bradley_ayanamsa(dt) if zodiac == "sidereal" else 0.0
        longs = {n: normalize_deg(L - (ay if zodiac=="sidereal" else 0.0)) for n, L in longs_trop.items()}

        # ASC/MC
        asc_trop, mc_trop = asc_mc(dt, lat, lon)
        asc = normalize_deg(asc_trop - (ay if zodiac=="sidereal" else 0.0))
        mc = normalize_deg(mc_trop - (ay if zodiac=="sidereal" else 0.0))

        # Houses (Equal)
        cusps = equal_house_cusps(asc, mode=house_mode)
        if house_clockwise:
            cusps = list(reversed(cusps))

        # Aspects
        aspects_found = find_aspects(longs, aspect_opts)
        aspect_orbs_dict = { spec['name']: aspect_opts.get(spec['key']+"_orb", spec['default_orb']) for spec in ASPECTS_DEF }

        # LST for header
        lst_val_deg = lst_deg(dt, lon)
        lst_hours = (lst_val_deg / 15.0) % 24.0
        lst_h = int(lst_hours); lst_m = int((lst_hours - lst_h) * 60); lst_s = int(round((((lst_hours - lst_h) * 60) - lst_m) * 60))
        if lst_s == 60:
            lst_s = 0; lst_m += 1
        if lst_m == 60:
            lst_m = 0; lst_h = (lst_h + 1) % 24
        lst_str = f"{lst_h:02d}:{lst_m:02d}:{lst_s:02d}"

        # Build tables
        table_rows = []
        for name in [n for n,_ in PLANETS if n != 'Earth']:
            if name not in longs:
                continue
            lam = longs[name]
            table_rows.append({"name": name, "lon": lam, "lon_fmt": format_longitude(lam)})

        houses_rows = [{"idx": i+1, "lon": c, "lon_fmt": format_longitude(c)} for i, c in enumerate(cusps)]

        data = {
            "person": person,
            "place": place,
            "dt_disp": dt_disp,
            "tz": tz_name,
            "utc_offset": utc_offset,
            "lst": lst_str,
            "lat": round(lat, 6),
            "lon": round(lon, 6),
            "elev": elev,
            "frame": frame,
            "zodiac": "Sidereal" if zodiac=="sidereal" else "Tropical",
            "ayanamsa": round(ay, 6) if zodiac=="sidereal" else 0.0,
            "asc_fmt": format_longitude(asc),
            "mc_fmt": format_longitude(mc),
            "houses": houses_rows,
            "table": table_rows,
            "aspects": aspects_found,
            "aspect_orbs": aspect_orbs_dict,
        }

        data_json = {
            "rotationDeg": normalize_deg(180.0 + asc),  # ASC at left
            "cusps": cusps,
            "rows": [{"name": r["name"], "lon": r["lon"]} for r in table_rows],
            "aspects": [{"p1": a['p1'], "p2": a['p2'], "type": a['type'], "delta": a['delta']} for a in aspects_found],
        }

        # Echo defaults/aspects for form re-render
        aspects_ui = [
            {"key": spec["key"], "name": spec["name"], "orb": aspect_opts.get(spec["key"]+"_orb", spec["default_orb"]), "on": aspect_opts.get(spec["key"]+"_on", True)}
            for spec in ASPECTS_DEF
        ]

        # Keep the user-entered birth time in the form
        default_dt = local_dt.strftime("%Y-%m-%dT%H:%M")

        return render_template_string(
            LAYOUT,
            default_dt=default_dt,
            default_tz=tz_name,
            data=data,
            data_json=data_json,
            aspects=aspects_ui,
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 400


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=10000)

