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
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>New Astrology Emerging — Switchboard (Local)</title>
  <style>
    :root { --bg:#0b0f14; --card:#121821; --ink:#eaf2ff; --muted:#9db2cf; --accent:#7cc0ff; }
    html,body { background:var(--bg); color:var(--ink); font-family: ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto; }
    .wrap { max-width: 1200px; margin: 24px auto; padding: 0 16px; }
    .card { background:var(--card); border-radius: 16px; padding: 18px; box-shadow: 0 10px 30px rgba(0,0,0,.25); }
    h1 { font-weight:700; letter-spacing:.2px; }
    label { display:block; margin:10px 0 6px; color:var(--muted); font-size:14px; }
    input, select { background:#0e1520; color:var(--ink); border:1px solid #1f2a38; border-radius:10px; padding:10px 12px; width:100%; }
    .row { display:grid; grid-template-columns: repeat(3, 1fr); gap: 14px; }
    .row2 { display:grid; grid-template-columns: repeat(2, 1fr); gap: 14px; }
    .row3 { display:grid; grid-template-columns: repeat(3, 1fr); gap: 14px; }
    .actions { display:flex; gap:10px; margin-top:14px; }
    button { background:var(--accent); color:#001; border:0; padding:12px 16px; border-radius:12px; font-weight:700; cursor:pointer; }
    .grid { display:grid; grid-template-columns: 520px 1fr; gap:16px; }
    canvas { background:#0c131e; border:1px solid #1f2a38; border-radius:16px; width:100%; height:auto; }
    table { width:100%; border-collapse: collapse; }
    th, td { text-align:left; padding:8px 6px; border-bottom:1px dashed #253246; }
    .muted { color:var(--muted); }
    .pill { display:inline-block; background:#0f1a28; border:1px solid #1e2a3a; padding:4px 8px; border-radius:999px; font-size:12px; }
    .section-title{ margin-top:10px; font-weight:700; }
    .minihead{ display:flex; gap:16px; align-items:center; font-size:14px; }
    .minihead b{ font-size:16px; }
    .dot::before{ content:"•"; margin:0 8px; color:#5a6e8a; }
  @media print {
    body { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .card, .wrap { box-shadow: none; }
    form, .actions, details, .report-text { display: none !important; }
    .minihead { display:none !important; }
    .grid { grid-template-columns: 1fr !important; }
    canvas { width: 7.5in !important; height: 7.5in !important; }
    table { font-size: 12px; }
    th, td { padding: 4px 6px; }
    .card { padding: 0; background: none; border: 0; }
  }
  /* page size for print */
  @page { size: letter; margin: 0.5in; }
    .card, .wrap { box-shadow: none; }
    form, .actions, details { display: none; }
    .grid { grid-template-columns: 1fr; }
    .report-chart { break-after: page; page-break-after: always; }
    .report-text { break-before: page; page-break-before: always; }
    canvas { width: 9in; height: 9in; }
  }
    .card, .wrap { box-shadow: none; }
    form, .actions, details { display: none; }
    .grid { grid-template-columns: 1fr; }
    canvas { width: 9in; height: 9in; }
  }
  @page { size: letter; margin: 0.5in; }
  .tooltip{ position:fixed; z-index:1000; background:#0e1520; color:#eaf2ff; border:1px solid #1f2a38; border-radius:8px; padding:6px 8px; font-size:13px; pointer-events:none; box-shadow:0 6px 20px rgba(0,0,0,.35); }
  </style>
</head>
<body>
  <div class="wrap">
    <h1>New Astrology Emerging — Switchboard (Local)</h1>
    <p class="muted">Birth charts with equal houses (Asc middle/cusp), custom aspect orbs, and a clean printable report.</p>

    <div class="card">
      <form method="GET" action="{{ url_for('chart') }}">
        <label>Name (optional)</label>
        <input name="person" placeholder="Full name" value="{{ data['person'] if data else '' }}">
        <div class="row2">
          <div>
            <label>Birth date & time</label>
            <input type="datetime-local" name="dt" value="{{ default_dt }}" required>
            <label style="margin-top:8px">Timezone</label>
            <select name="tz">
              {% if data %}
              <option value="{{ data['tz'] }}" selected>{{ data['tz'] }}{% if request.args.get('tz','auto')=='auto' %} (auto){% endif %}</option>
              {% endif %}
              <option value="auto" {% if (data and data['tz']=='auto') or (not data and default_tz=='auto') %}selected{% endif %}>Auto (from birthplace)</option>
              <option value="UTC" {% if (data and data['tz']=='UTC') or (not data and default_tz=='UTC') %}selected{% endif %}>UTC</option>
              <option value="America/Denver" {% if (data and data['tz']=='America/Denver') or (not data and default_tz=='America/Denver') %}selected{% endif %}>America/Denver (MT)</option>
              <option value="America/Los_Angeles" {% if (data and data['tz']=='America/Los_Angeles') or (not data and default_tz=='America/Los_Angeles') %}selected{% endif %}>America/Los_Angeles (PT)</option>
              <option value="America/New_York" {% if (data and data['tz']=='America/New_York') or (not data and default_tz=='America/New_York') %}selected{% endif %}>America/New_York (ET)</option>
              <option value="Europe/London" {% if (data and data['tz']=='Europe/London') or (not data and default_tz=='Europe/London') %}selected{% endif %}>Europe/London</option>
              <option value="Europe/Paris" {% if (data and data['tz']=='Europe/Paris') or (not data and default_tz=='Europe/Paris') %}selected{% endif %}>Europe/Paris</option>
              <option value="Asia/Kolkata" {% if (data and data['tz']=='Asia/Kolkata') or (not data and default_tz=='Asia/Kolkata') %}selected{% endif %}>Asia/Kolkata</option>
              <option value="Asia/Tokyo" {% if (data and data['tz']=='Asia/Tokyo') or (not data and default_tz=='Asia/Tokyo') %}selected{% endif %}>Asia/Tokyo</option>
              <option value="Australia/Sydney" {% if (data and data['tz']=='Australia/Sydney') or (not data and default_tz=='Australia/Sydney') %}selected{% endif %}>Australia/Sydney</option>
            </select>
          </div>
          <div>
            <label>Birthplace (city, region, country)</label>
            <input name="place" placeholder="e.g., Boulder, CO, USA" value="{{ data['place'] if data else '' }}">
            <details style="margin-top:8px"><summary class="muted">Advanced: Coordinates (optional)</summary>
              <div class="row">
                <input name="lat" placeholder="Latitude e.g., 40.014986" value="{{ data['lat'] if data else '' }}">
                <input name="lon" placeholder="Longitude e.g., -105.270546" value="{{ data['lon'] if data else '' }}">
                <input name="elev" placeholder="Elevation m (optional)" value="{{ data['elev'] if data else '' }}">
              </div>
              <p class="muted">If provided, coordinates override the birthplace lookup.</p>
            </details>
          </div>
        </div>
        <div class="row3">
          <div>
            <label>Frame</label>
            <select name="frame">
              <option value="geo" selected>Geocentric</option>
              <option value="helio">Heliocentric</option>
            </select>
          </div>
          <div>
            <label>Zodiac</label>
            <select name="zodiac">
              <option value="tropical">Tropical</option>
              <option value="sidereal" selected>Sidereal (Fagan/Bradley)</option>
            </select>
          </div>
          <div>
            <label>House system</label>
            <select name="house_mode">
              <option value="asc_middle" selected>Equal – Asc in middle of 1st</option>
              <option value="asc_cusp">Equal – Asc at cusp of 1st</option>
            </select>
          </div>
        </div>
        <div class="row3">
          <div>
            <label>Clockwise house numbering</label>
            <select name="house_clockwise">
              <option value="no" selected>No (standard)</option>
              <option value="yes">Yes (1→12 clockwise)</option>
            </select>
          </div>
          <div>
            <label>Aspect set</label>
            <div class="row" style="grid-template-columns: repeat(6,1fr);">
              {% for a in aspects %}
              <label style="display:flex;align-items:center;gap:6px;margin:0"><input type="checkbox" name="{{a.key}}_on" {% if a.on %}checked{% endif %}> {{a.name}}</label>
              {% endfor %}
            </div>
          </div>
          <div>
            <label>Orbs (°)</label>
            <div class="row" style="grid-template-columns: repeat(6,1fr);">
              {% for a in aspects %}
                <input type="number" step="0.1" min="0" max="15" name="{{a.key}}_orb" value="{{a.orb}}" placeholder="{{a.name}} orb">
              {% endfor %}
            </div>
          </div>
        </div>
        <div class="actions">
          <button type="submit">Compute Chart</button>
          <a class="pill" href="{{ url_for('about') }}">About & limits</a>
        </div>
      </form>
    </div>

    {% if data %}
    <div class="card minihead">
      <div><b>{{ data['person'] or '—' }}</b></div>
      <div class="dot"></div>
      <div>{{ data['place'] or '—' }}</div>
      <div class="dot"></div>
      <div>{{ data['dt_disp'] }} ({{ data['tz'] }}, UTC{{ data['utc_offset'] }})</div>
      <div class="dot"></div>
      <div>LST {{ data['lst'] }}</div>
      <div class="dot"></div>
      <div>Lat {{ data['lat'] }}°, Lon {{ data['lon'] }}°</div>
    </div>
    <div class="report-chart">
    <div class="grid" style="margin-top:16px;">
      <div class="card">
        <canvas id="wheel" width="860" height="860"></canvas>
        <div id="chartTip" class="tooltip" style="display:none"></div>
      </div>
      <div class="card">
        <h3>Positions ({{ data['frame'].upper() }} · {{ data['zodiac'].title() }})</h3>
        <div class="row2">
          <div>
            <table>
              <thead><tr><th>Body</th><th>Longitude</th></tr></thead>
              <tbody>
                {% for row in data['table'] %}
                  <tr><td>{{ row.name }}</td><td>{{ row.lon_fmt }}</td></tr>
                {% endfor %}
              </tbody>
            </table>
            <p class="muted">ASC: {{ data['asc_fmt'] }} · MC: {{ data['mc_fmt'] }}</p>
            <p class="muted">Place: {{ data['place'] or '—' }} ({{ data['lat'] }}, {{ data['lon'] }}) · TZ: {{ data['tz'] }}</p>
            <p class="muted">Ayanamsa = {{ data['ayanamsa'] }}° (Fagan/Bradley approx). Time: {{ data['dt_disp'] }} (UTC{{ data['utc_offset'] }})</p>

            <details style="margin-top:8px"><summary class="muted">Notes</summary>
              <ul class="muted">
                <li>Birth time accuracy matters — ASC/MC and houses shift quickly.</li>
                <li>“Auto” timezone comes from birthplace; override it if it looks wrong.</li>
                <li>If geocoding can’t find the city, enter coordinates under Advanced.</li>
                <li>Sidereal ayanamsa uses Swiss Ephemeris when installed (fallback to polynomial).</li>
              </ul>
            </details>
          </div>
          <div>
            <div class="section-title">Houses (Equal)</div>
            <table>
              <thead><tr><th>#</th><th>Cusp</th></tr></thead>
              <tbody>
              {% for h in data['houses'] %}
                <tr><td>{{ h.idx }}</td><td>{{ h.lon_fmt }}</td></tr>
              {% endfor %}
              </tbody>
            </table>
          </div>
        </div>

        <div class="section-title">Aspects</div>
        <table>
          <thead><tr><th>Aspect</th><th>Between</th><th>Exact</th><th>Orb (°)</th><th>Δ</th></tr></thead>
          <tbody>
            {% for a in data['aspects'] %}
            <tr>
              <td>{{ a.type }}</td>
              <td>{{ a.p1 }} – {{ a.p2 }}</td>
              <td>{{ a.angle }}°</td>
              <td>{{ data['aspect_orbs'][a.type] }}</td>
              <td>{{ a.delta }}°</td>
            </tr>
            {% endfor %}
            {% if data['aspects']|length == 0 %}
              <tr><td colspan="5" class="muted">No aspects within chosen orbs.</td></tr>
            {% endif %}
          </tbody>
        </table>

        <div class="section-title">Orbs Used</div>
        <table>
          <thead><tr><th>Aspect</th><th>Orb (°)</th></tr></thead>
          <tbody>
            {% for a in aspects %}
              <tr><td>{{ a.name }}</td><td>{{ a.orb }}</td></tr>
            {% endfor %}
          </tbody>
        </table>

        

        <div class="actions">
          <button onclick="window.print()">Print / Save PDF</button>
        </div>
      </div>
    </div>
  </div>

  <div class="card report-text">
    <h3>Text listing</h3>
    <pre class="muted">Planets:
{% for row in data['table'] %} - {{ row.name }}: {{ row.lon_fmt }}
{% endfor %}Aspects:
{% if data['aspects']|length == 0 %} - (none within chosen orbs)
{% else %}
{% for a in data['aspects'] %} - {{ a.type }}: {{ a.p1 }} – {{ a.p2 }} (Δ {{ a.delta }}°)
{% endfor %}
{% endif %}
    </pre>
  </div>

    <script>
  // Chart rendering with header + clearer planet styling
  const table = {{ data_json | tojson }};
  const header = {{ {'person': data['person'], 'place': data['place'], 'dt': data['dt_disp'], 'tz': data['tz'], 'offset': data['utc_offset'], 'lat': data['lat'], 'lon': data['lon'], 'lst': data['lst']} | tojson }};

  const canvas = document.getElementById('wheel');
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height; const cx=W/2, cy=H/2;
  const R = Math.min(W,H)*0.45;

  const signGlyph = ['♈','♉','♊','♋','♌','♍','♎','♏','♐','♑','♒','♓'];
  const planetGlyph = { Sun:'☉', Moon:'☾', Mercury:'☿', Venus:'♀', Mars:'♂', Jupiter:'♃', Saturn:'♄', Uranus:'♅', Neptune:'♆', Pluto:'♇' };
  const planetColor = { Sun:'#FFD166', Moon:'#B0BEC5', Mercury:'#64B5F6', Venus:'#F48FB1', Mars:'#EF5350', Jupiter:'#FFB74D', Saturn:'#D7CCC8', Uranus:'#4DD0E1', Neptune:'#64B5F6', Pluto:'#BA68C8' };

  function deg2rad(d){ return d*Math.PI/180; }
  function drawCircle(r, w=2, stroke='#2a3a52'){ ctx.beginPath(); ctx.lineWidth=w; ctx.arc(cx,cy,r,0,Math.PI*2); ctx.strokeStyle=stroke; ctx.stroke(); }
  function drawTick(angleDeg, r1, r2, lw=1, stroke='#1f2a38'){
    const a = deg2rad(angleDeg), ca=Math.cos(a), sa=Math.sin(a);
    ctx.beginPath(); ctx.moveTo(cx+ca*r1, cy+sa*r1); ctx.lineTo(cx+ca*r2, cy+sa*r2);
    ctx.lineWidth = lw; ctx.strokeStyle = stroke; ctx.stroke();
  }
  function drawTextOnRing(txt, angleDeg, radius, font='18px ui-sans-serif', fill='#eaf2ff'){
    const a = deg2rad(angleDeg);
    const x = cx + Math.cos(a)*radius, y = cy + Math.sin(a)*radius;
    ctx.save(); ctx.translate(x,y); ctx.rotate(a + Math.PI/2);
    ctx.fillStyle = fill; ctx.font = font; ctx.textAlign='center'; ctx.textBaseline='middle';
    ctx.shadowColor='rgba(0,0,0,.45)'; ctx.shadowBlur=3; ctx.fillText(txt, 0, 0); ctx.restore();
  }
  function pad2(n){ return (n<10?'0':'')+n; }
  function degMin(d){ let a=((d%360)+360)%360; const D=Math.floor(a%30); const M=Math.floor((a-Math.floor(a))*60); return `${pad2(D)}°${pad2(M)}′`; }
  function roundRect(x, y, w, h, r){
    const rr = Math.min(r, Math.min(w,h)/2);
    ctx.beginPath(); ctx.moveTo(x+rr, y); ctx.arcTo(x+w, y, x+w, y+h, rr); ctx.arcTo(x+w, y+h, x, y+h, rr); ctx.arcTo(x, y+h, x, y, rr); ctx.arcTo(x, y, x+w, y, rr); ctx.closePath();
  }

  // Clear, base rings
  ctx.clearRect(0,0,W,H);
  drawCircle(R,3); drawCircle(R*0.92,1); drawCircle(R*0.78,1); drawCircle(R*0.64,1);

  const rotation = table.rotationDeg; // ASC at left

  // Header inside the canvas
  (function drawHeader(){
    const parts = [];
    if (header.person && header.person.trim() !== '') parts.push(header.person.trim());
    if (header.place && header.place.trim() !== '') parts.push(header.place.trim());
    parts.push(`${header.dt} (${header.tz}, UTC${header.offset})`);
    parts.push(`LST ${header.lst}`);
    parts.push(`Lat ${header.lat}°, Lon ${header.lon}°`);
    const title = parts.join('  •  ');
    ctx.save();
    ctx.fillStyle='#eaf2ff'; ctx.font='20px ui-sans-serif'; ctx.textAlign='center'; ctx.textBaseline='alphabetic';
    ctx.shadowColor='rgba(0,0,0,.6)'; ctx.shadowBlur=4; ctx.lineWidth=3;
    ctx.fillText(title, cx, cy - R - 24);
    ctx.restore();
  })();

  // 30° majors + 5° minors
  for(let d=0; d<360; d+=5){
    const major = (d%30===0);
    drawTick(-(d)+rotation, R, R*(major?0.94:0.97), major?2:1, major?'#32445f':'#1f2a38');
  }

  // Sign glyphs
  for (let i=0;i<12;i++){
    const mid=i*30+15; drawTextOnRing(signGlyph[i], -(mid)+rotation, R*1.02, '30px ui-sans-serif');
  }

  // Houses: rays + numbers
  table.cusps.forEach((cusp,i)=>{
    drawTick(-(cusp)+rotation, 0, R, 1.75, '#203045');
    drawTextOnRing(String(i+1), -(cusp+15)+rotation, R*0.88, '18px ui-sans-serif', '#cfe0ff');
  });

  // Precompute planet positions first (so we can draw aspects underneath)
  const positions = {};
  table.rows.forEach(row=>{
    const lon = row.lon;
    const a = deg2rad(-(lon)+rotation);
    const pr = R*0.80;
    const x = cx + Math.cos(a)*pr, y = cy + Math.sin(a)*pr;
    positions[row.name] = { x, y, lon };
  });

  // Aspect lines UNDER the planets for clarity
  const aspectStyle = { Conjunction:'#eaf2ff', Opposition:'#ff9aa2', Trine:'#9ae59a', Square:'#ffd49a', Sextile:'#9ac7ef', Quincunx:'#c9a4ef' };
  table.aspects.forEach(a=>{
    const p1 = positions[a.p1], p2 = positions[a.p2];
    if(!p1||!p2) return;
    ctx.beginPath(); ctx.moveTo(p1.x,p1.y); ctx.lineTo(p2.x,p2.y);
    ctx.strokeStyle = aspectStyle[a.type] || '#a4b3c6';
    ctx.lineWidth = 2; ctx.globalAlpha = 0.9; ctx.stroke(); ctx.globalAlpha = 1;
  });

  // Planets — GLYPH ONLY (larger, high-contrast)
  table.rows.forEach(row=>{
    const p = positions[row.name];
    const col = planetColor[row.name] || '#7cc0ff';
    // soft halo + colored disc
    ctx.save(); ctx.shadowColor = col; ctx.shadowBlur = 12;
    ctx.beginPath(); ctx.arc(p.x, p.y, 18, 0, Math.PI*2); ctx.fillStyle = col; ctx.fill(); ctx.restore();
    // outline for separation
    ctx.beginPath(); ctx.arc(p.x, p.y, 18, 0, Math.PI*2); ctx.lineWidth = 2; ctx.strokeStyle = '#0b0f14'; ctx.stroke();
    // glyph centered, with dark outline + white fill
    const g = planetGlyph[row.name] || row.name[0];
    ctx.font = '26px ui-sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
    ctx.strokeStyle = '#0b0f14'; ctx.lineWidth = 3; ctx.strokeText(g, p.x, p.y);
    ctx.fillStyle = '#ffffff'; ctx.fillText(g, p.x, p.y);
  });

  // === Hover tooltips for planets and zodiac ring ===
  const tipEl = document.getElementById('chartTip');
  const signNames = ['Aries','Taurus','Gemini','Cancer','Leo','Virgo','Libra','Scorpio','Sagittarius','Capricorn','Aquarius','Pisces'];
  function norm360(a){ return ((a%360)+360)%360; }
  function showTip(html, clientX, clientY){
    tipEl.innerHTML = html; tipEl.style.display='block';
    const pad=14; tipEl.style.left = (clientX + pad) + 'px'; tipEl.style.top = (clientY + pad) + 'px';
  }
  function hideTip(){ tipEl.style.display='none'; }

  canvas.addEventListener('mouseleave', hideTip);
  canvas.addEventListener('mousemove', (e)=>{
    const rect = canvas.getBoundingClientRect();
    // account for CSS scaling vs canvas pixels
    const scaleX = canvas.width / rect.width, scaleY = canvas.height / rect.height;
    const mx = (e.clientX - rect.left) * scaleX;
    const my = (e.clientY - rect.top) * scaleY;

    const dx = mx - cx, dy = my - cy; const r = Math.hypot(dx,dy);
    const ang = Math.atan2(dy, dx); // radians
    const lon = norm360(-ang*180/Math.PI + table.rotationDeg);

    // 1) Planet hit test (within 22px of planet center)
    let closest = null, cd = 1e9;
    table.rows.forEach(row=>{
      const p = positions[row.name]; if(!p) return;
      const d = Math.hypot(mx - p.x, my - p.y);
      if (d < cd){ cd = d; closest = {row, p}; }
    });
    if (closest && cd <= 22){
      const name = closest.row.name;
      const sIdx = Math.floor(norm360(closest.p.lon)/30);
      const html = `<b>${planetGlyph[name]||name}</b> ${name}<br>${signGlyph[sIdx]} ${signNames[sIdx]} · ${degMin(closest.p.lon)}`;
      showTip(html, e.clientX, e.clientY); return;
    }

    // 2) Zodiac ring hover (near outer ring)
    if (Math.abs(r - R) < 24){
      const sIdx = Math.floor(lon/30);
      const html = `${signGlyph[sIdx]} <b>${signNames[sIdx]}</b><br>${degMin(lon)}`;
      showTip(html, e.clientX, e.clientY); return;
    }

    hideTip();
  });
</script>
    {% endif %}
  </div>
</body>
</html>
"""

ABOUT = """
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

@app.route("/")
def index():
    # Default dropdown selection uses 'auto' (detect from birthplace)
    default_tz = 'auto'
    # Prefill the datetime-local field with a valid string (any zone works just to seed the control)
    now_local = datetime.now(pytz.timezone('America/Denver')).replace(second=0, microsecond=0)
    default_dt = now_local.strftime("%Y-%m-%dT%H:%M")

    # aspect controls default state
    aspects = []
    for spec in ASPECTS_DEF:
        aspects.append({"key": spec['key'], "name": spec['name'], "orb": spec['default_orb'], "on": True})

    return render_template_string(LAYOUT, default_dt=default_dt, default_tz=default_tz, data=None, aspects=aspects)

@app.route("/about")
def about():
    return ABOUT

@app.route("/chart")
def chart():
    try:
        # Core fields
        dt_str = request.args.get("dt")
        frame = request.args.get("frame", "geo")
        zodiac = request.args.get("zodiac", "sidereal")
        house_mode = request.args.get("house_mode", "asc_middle")
        house_clockwise = request.args.get("house_clockwise", "no") == "yes"
        person = (request.args.get("person") or "").strip()

        # Aspect form inputs
        aspect_opts = {}
        for spec in ASPECTS_DEF:
            aspect_opts[spec['key']+"_on"] = (request.args.get(spec['key']+"_on") is not None)
            orb_val = request.args.get(spec['key']+"_orb")
            aspect_opts[spec['key']+"_orb"] = float(orb_val) if orb_val not in (None, '') else spec['default_orb']

        # Birthplace or coordinates
        place = (request.args.get("place") or "").strip()
        lat_str = (request.args.get("lat") or "").strip()
        lon_str = (request.args.get("lon") or "").strip()
        elev_str = (request.args.get("elev") or "").strip()

        lat = lon = None
        elev = float(elev_str) if elev_str else 0.0

        if lat_str and lon_str:
            lat = float(lat_str); lon = float(lon_str)
            resolved_place = place or "(coords provided)"
        elif place:
            try:
                loc = _geocoder.geocode(place, addressdetails=True, language='en', timeout=10)
                if not loc:
                    return jsonify({'error': f"Couldn't find that birthplace: '{place}'. Try a broader query or enter coordinates."}), 400
                lat = float(loc.latitude); lon = float(loc.longitude)
                resolved_place = getattr(loc, 'address', getattr(loc, 'display_name', place))
            except Exception as ge:
                return jsonify({'error': f"Geocoding failed: {ge}"}), 400
        else:
            return jsonify({'error': 'Please enter a birthplace or coordinates.'}), 400

        # Timezone handling (auto from coords, or manual)
        tz_name = request.args.get("tz", "auto")
        if tz_name == 'auto':
            try:
                guess = _tzf.timezone_at(lng=lon, lat=lat)
            except Exception:
                guess = None
            tz_name = guess or 'UTC'
        try:
            tz = pytz.timezone(tz_name)
        except Exception:
            tz = pytz.UTC
            tz_name = 'UTC'

        # Localize birth time and convert to UTC
        local_naive = datetime.strptime(dt_str, "%Y-%m-%dT%H:%M")
        local_dt = tz.localize(local_naive)
        dt = local_dt.astimezone(pytz.UTC)

        # For display: UTC offset like +HH:MM
        offset_sec = int(local_dt.utcoffset().total_seconds())
        sign = '+' if offset_sec >= 0 else '-'
        hh = abs(offset_sec)//3600
        mm = (abs(offset_sec)%3600)//60
        utc_offset = f"{sign}{hh:02d}:{mm:02d}"

        helio = (frame == 'helio')
        longs_trop = planetary_longitudes(dt, lat, lon, elev, helio=helio)

        if zodiac == 'sidereal':
            longs = to_sidereal(longs_trop, dt)
            ay = round(fagan_bradley_ayanamsa(dt), 6)
        else:
            longs = longs_trop
            ay = 0.0

        # ASC / MC and houses (computed in same zodiac as positions)
        asc_trop, mc_trop = asc_mc(dt, lat, lon)
        if zodiac == 'sidereal':
            asc = normalize_deg(asc_trop - ay)
            mc = normalize_deg(mc_trop - ay)
        else:
            asc, mc = asc_trop, mc_trop

        cusps = equal_house_cusps(asc, mode=house_mode)

        # Local Sidereal Time (at birthplace longitude)
        lst_val_deg = lst_deg(dt, lon)
        lst_hours = (lst_val_deg / 15.0) % 24.0
        lst_h = int(lst_hours)
        lst_m = int((lst_hours - lst_h) * 60)
        lst_s = int(round((((lst_hours - lst_h) * 60) - lst_m) * 60))
        if lst_s == 60:
            lst_s = 0
            lst_m += 1
        if lst_m == 60:
            lst_m = 0
            lst_h = (lst_h + 1) % 24
        lst_str = f"{lst_h:02d}:{lst_m:02d}:{lst_s:02d}"

        rows = []
        for name in [k for k,_ in PLANETS if (k != 'Earth' or helio)]:
            lonv = longs.get(name)
            rows.append({
                'name': name,
                'lon': round(lonv, 6),
                'lon_str': format_longitude(lonv),
            })

        # Aspects
        aspects_list = find_aspects(longs, aspect_opts)

        aspect_orbs_dict = { spec['name']: aspect_opts.get(spec['key']+"_orb", spec['default_orb']) for spec in ASPECTS_DEF }

        # House label order (numbering)
        labels = [f"{i}" for i in range(1,13)]  # visual orientation is handled by rotation

        data = {
            'frame': frame,
            'zodiac': zodiac,
            'ayanamsa': round(ay, 6),
            'dt_disp': local_dt.strftime('%Y-%m-%d %H:%M'),
            'utc_offset': utc_offset,
            'tz': tz_name,
            'place': resolved_place,
            'lat': round(lat, 6),
            'lon': round(lon, 6),
            'elev': elev,
            'lst': lst_str,
            'table': [{'name': r['name'], 'lon_fmt': r['lon_str']} for r in rows],
            'asc_fmt': format_longitude(asc),
            'mc_fmt': format_longitude(mc),
            'houses': [{'idx': (i+1), 'lon_fmt': format_longitude(l)} for i, l in enumerate(cusps)],
            'aspects': aspects_list,
            'person': person,
            'aspect_orbs': aspect_orbs_dict,
        }

        data_json = {
            'rows': rows,
            'asc': asc,
            'mc': mc,
            'rotationDeg': asc + 180.0,  # put ASC at left
            'cusps': cusps,
            'houseClockwise': house_clockwise,
            'houseLabels': labels,
            'aspects': aspects_list,
        }

        # keep the user-entered birth time in the form
        default_dt = local_dt.strftime("%Y-%m-%dT%H:%M")

        aspects_ui = []
        for spec in ASPECTS_DEF:
            aspects_ui.append({
                'key': spec['key'], 'name': spec['name'],
                'orb': aspect_opts.get(spec['key']+"_orb", spec['default_orb']),
                'on': aspect_opts.get(spec['key']+"_on", True)
            })

        return render_template_string(
            LAYOUT,
            default_dt=default_dt,
            default_tz=tz_name,
            data=data,
            data_json=data_json,
            aspects=aspects_ui
        )

    except Exception as e:
        return jsonify({ 'error': str(e) }), 400

if __name__ == "__main__":
    app.run(debug=True)
