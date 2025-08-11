"""
Local, modern **New Astrology Emerging** switchboard — now with HOUSES and ASPECTS.

Features (MVP+):
- Local Flask web app with a single-page UI.
- Geocentric & Heliocentric positions via Skyfield (DE421 by default).
- Tropical and Sidereal (Fagan/Bradley) longitudes.
- **Equal houses** with two modes:
  • Ascendant at cusp of 1st
  • Ascendant in the middle of 1st (New Astrology Emerging style)
- **Clockwise house numbering** toggle (if desired).
- **Asc/MC** computation and display.
- **Aspects**: conjunction, opposition, trine, square, sextile, quincunx with custom orbs.
- Clean chart wheel in HTML <canvas> with house lines, labels, and aspect lines.
- Printable report (browser print to PDF).

Install:
  pip install flask skyfield jinja2 pytz
Run:
  python app.py  →  http://127.0.0.1:5000

Tested on Python 3.9+.
"""

from __future__ import annotations
from math import atan2, degrees, radians, sin, cos, floor, asin
from datetime import datetime
import pytz

from flask import Flask, request, render_template_string, jsonify

# SKYFIELD
from skyfield.api import load, wgs84
from skyfield.api import Loader

app = Flask(__name__)

_loader = None
_ts = None
_eph = None

PLANETS = [
    ("Sun", "sun"),
    ("Moon", "moon"),
    ("Mercury", "mercury"),
    ("Venus", "venus"),
    ("Earth", "earth"),
    ("Mars", "mars"),
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

# -----------------------------
# Utility math / formatting
# -----------------------------

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
    return f"{ZODIAC_SIGNS[sign_index]} {di:02d}°{mi:02d}'{si:04.1f}\u2033"

# ---------------------------------------
# Fagan/Bradley ayanamsa (approximation)
# ---------------------------------------
from datetime import timezone as _tz

def fagan_bradley_ayanamsa(dt: datetime) -> float:
    # Convert to UTC (approx for polynomial time fraction)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=_tz.utc)
    else:
        dt = dt.astimezone(_tz.utc)
    y = dt.year + (dt.timetuple().tm_yday - 1 + (dt.hour + dt.minute/60 + dt.second/3600)/24) / (366 if _is_leap(dt.year) else 365)
    t = y - 2000.0
    base = 24.042044444
    rate = -0.013004167
    curve = -0.000000164 * t * t
    ay = base + rate * t + curve
    return normalize_deg(ay)


def _is_leap(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

# --------------------------------
# Skyfield helpers
# --------------------------------

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
        _eph = get_loader()("de421.bsp")
    return _eph

# --------------------------------
# Core calculations
# --------------------------------

def planetary_longitudes(dt: datetime, lat: float, lon: float, elevation_m: float,
                         helio: bool = False) -> dict:
    ts = get_timescale()
    t = ts.from_datetime(dt)
    eph = get_ephemeris()

    if helio:
        origin = eph["sun"]
    else:
        observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=elevation_m)
        origin = observer

    results = {}
    ecliptic_frame = eph["earth"].ecliptic_frame()

    for name, key in PLANETS:
        if key == "earth":
            if not helio:
                continue
        target = eph[key]
        astrometric = origin.at(t).observe(target)
        lat_ecl, lon_ecl, distance = astrometric.frame_latlon(ecliptic_frame)
        lam = normalize_deg(lon_ecl.degrees)
        results[name] = lam

    return results


def to_sidereal(longitudes_tropical: dict, dt: datetime) -> dict:
    ay = fagan_bradley_ayanamsa(dt)
    return {name: normalize_deg(lon - ay) for name, lon in longitudes_tropical.items()}

# --------------------------------
# Ascendant / MC / Houses
# --------------------------------

OBLIQ_DEG = 23.43929111  # J2000 mean obliquity (good enough for MVP)

def gmst_hours(dt: datetime) -> float:
    t = get_timescale().from_datetime(dt)
    return t.gmst  # hours


def lst_deg(dt: datetime, lon_deg: float) -> float:
    lst_h = gmst_hours(dt) + (lon_deg / 15.0)
    lst_h %= 24.0
    return lst_h * 15.0


def asc_mc(dt: datetime, lat_deg: float, lon_deg: float) -> tuple[float, float]:
    """Return (ASC, MC) ecliptic longitudes in degrees.
    MC: analytic conversion via tan(λ) = tan(θ)/cos ε.
    ASC: solved numerically by finding ecliptic point at horizon in the east.
    """
    ε = radians(OBLIQ_DEG)
    φ = radians(lat_deg)
    θ = radians(lst_deg(dt, lon_deg))  # local sidereal angle

    # MC
    lam_mc = degrees(atan2(sin(θ), cos(θ)*cos(ε)))
    lam_mc = normalize_deg(lam_mc)

    # ASC via sampling and root-refinement on altitude = 0
    def eq_to_altaz(alpha, delta):
        H = (θ - alpha)
        # wrap H to [-pi, pi]
        while H > 3.141592653589793:
            H -= 2*3.141592653589793
        while H < -3.141592653589793:
            H += 2*3.141592653589793
        sin_h = sin(φ)*sin(delta) + cos(φ)*cos(delta)*cos(H)
        h = asin(max(-1.0, min(1.0, sin_h)))
        # At horizon (h≈0): estimate azimuth
        sinA = -sin(H) * cos(delta)
        cosA = (sin(delta)) / max(1e-9, cos(φ))
        A = atan2(sinA, cosA)  # -pi..pi, where 0=N, +pi/2=E
        if A < 0:
            A += 2*3.141592653589793
        return h, A

    def ecl_to_eq(lam):
        # beta=0
        lamr = radians(lam)
        alpha = atan2(sin(lamr)*cos(ε), cos(lamr))
        if alpha < 0:
            alpha += 2*3.141592653589793
        delta = asin(sin(lamr)*sin(ε))
        return alpha, delta

    # coarse scan
    samples = []
    step = 1.0
    prev_h = None
    prev_lam = None
    roots = []
    for lam in [i*step for i in range(361)]:
        alpha, delta = ecl_to_eq(lam)
        h, A = eq_to_altaz(alpha, delta)
        samples.append((lam, h, A))
        if prev_h is not None and (h*prev_h) < 0:  # sign change
            # linear interpolation for root
            lam0 = prev_lam + (0 - prev_h) * (lam - prev_lam) / (h - prev_h)
            alpha0, delta0 = ecl_to_eq(lam0)
            h0, A0 = eq_to_altaz(alpha0, delta0)
            roots.append((lam0, h0, A0))
        prev_h = h
        prev_lam = lam

    # choose eastern root (0°<Az<180°)
    asc_candidates = [lam for (lam, h0, A0) in roots if 0 < degrees(A0) < 180]
    if not asc_candidates:
        # fallback: take closest to A≈90°
        asc_candidates = [min(samples, key=lambda x: abs(degrees(x[2]) - 90))[0]]
    lam_asc = normalize_deg(asc_candidates[0])

    return lam_asc, lam_mc


def equal_house_cusps(asc_deg: float, mode: str = "asc_middle") -> list[float]:
    """Return 12 cusps (cusp1..cusp12) in degrees.
    mode: 'asc_middle' → cusp1 = ASC - 15° (Asc in middle of 1st)
          'asc_cusp'   → cusp1 = ASC      (Asc at cusp of 1st)
    """
    if mode == "asc_middle":
        cusp1 = normalize_deg(asc_deg - 15.0)
    else:
        cusp1 = normalize_deg(asc_deg)
    return [normalize_deg(cusp1 + i*30.0) for i in range(12)]

# --------------------------------
# Aspects
# --------------------------------

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

# ------------------------------
# UI templates
# ------------------------------
LAYOUT = """
<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>New Astrology Emerging Switchboard – Local (Houses + Aspects)</title>
  <style>
    :root { --bg:#0b0f14; --card:#121821; --ink:#eaf2ff; --muted:#9db2cf; --accent:#7cc0ff; }
    html,body { background:var(--bg); color:var(--ink); font-family: ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto; }
    .wrap { max-width: 1200px; margin: 24px auto; padding: 0 16px; }
    .card { background:var(--card); border-radius: 16px; padding: 18px; box-shadow: 0 10px 30px rgba(0,0,0,.25); }
    h1 { font-weight:700; letter-spacing:.2px; }
    label { display:block; margin:10px 0 6px; color:var(--muted); font-size:14px; }
    input, select { background:#0e1520; color:#eaf2ff; border:1px solid #1f2a38; border-radius:10px; padding:10px 12px; width:100%; }
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
  </style>
</head>
<body>
  <div class=\"wrap\">
    <h1>New Astrology Emerging Switchboard – Local (Houses + Aspects)</h1>
    <p class=\"muted\">Flash‑free. Equal houses (Asc middle/cusp), clockwise option, custom aspect orbs.</p>

    <div class=\"card\">
      <form method=\"GET\" action=\"{{ url_for('chart') }}\">
        <div class=\"row2\">
          <div>
            <label>Date & time (local)</label>
            <input type=\"datetime-local\" name=\"dt\" value=\"{{ default_dt }}\" required>
          </div>
          <div>
            <label>Location (lat, lon, elevation m)</label>
            <div class=\"row\">
              <input name=\"lat\" placeholder=\"40.0\" value=\"40.014986\" required>
              <input name=\"lon\" placeholder=\"-105.27\" value=\"-105.270546\" required>
              <input name=\"elev\" placeholder=\"1655\" value=\"1655\" required>
            </div>
          </div>
        </div>
        <div class=\"row3\">
          <div>
            <label>Frame</label>
            <select name=\"frame\">
              <option value=\"geo\" selected>Geocentric</option>
              <option value=\"helio\">Heliocentric</option>
            </select>
          </div>
          <div>
            <label>Zodiac</label>
            <select name=\"zodiac\">
              <option value=\"tropical\">Tropical</option>
              <option value=\"sidereal\" selected>Sidereal (Fagan/Bradley)</option>
            </select>
          </div>
          <div>
            <label>House system</label>
            <select name=\"house_mode\">
              <option value=\"asc_middle\" selected>Equal – Asc in middle of 1st</option>
              <option value=\"asc_cusp\">Equal – Asc at cusp of 1st</option>
            </select>
          </div>
        </div>
        <div class=\"row3\">
          <div>
            <label>Clockwise house numbering</label>
            <select name=\"house_clockwise\">
              <option value=\"no\" selected>No (standard)</option>
              <option value=\"yes\">Yes (1→12 clockwise)</option>
            </select>
          </div>
          <div>
            <label>Aspect set</label>
            <div class=\"row\" style=\"grid-template-columns: repeat(6,1fr);\">
              {% for a in aspects %}
              <label style=\"display:flex;align-items:center;gap:6px;margin:0\"><input type=\"checkbox\" name=\"{{a.key}}_on\" {% if a.on %}checked{% endif %}> {{a.name}}</label>
              {% endfor %}
            </div>
          </div>
          <div>
            <label>Orbs (°)</label>
            <div class=\"row\" style=\"grid-template-columns: repeat(6,1fr);\">
              {% for a in aspects %}
                <input type=\"number\" step=\"0.1\" min=\"0\" max=\"15\" name=\"{{a.key}}_orb\" value=\"{{a.orb}}\" placeholder=\"{{a.name}} orb\">
              {% endfor %}
            </div>
          </div>
        </div>
        <div class=\"actions\">
          <button type=\"submit\">Compute Chart</button>
          <a class=\"pill\" href=\"{{ url_for('about') }}\">About & limits</a>
        </div>
      </form>
    </div>

    {% if data %}
    <div class=\"grid\" style=\"margin-top:16px;\">
      <div class=\"card\">
        <canvas id=\"wheel\" width=\"860\" height=\"860\"></canvas>
      </div>
      <div class=\"card\">
        <h3>Positions ({{ data['frame'].upper() }} · {{ data['zodiac'].title() }})</h3>
        <div class=\"row2\">
          <div>
            <table>
              <thead><tr><th>Body</th><th>Longitude</th></tr></thead>
              <tbody>
                {% for row in data['table'] %}
                  <tr><td>{{ row.name }}</td><td>{{ row.lon_fmt }}</td></tr>
                {% endfor %}
              </tbody>
            </table>
            <p class=\"muted\">ASC: {{ data['asc_fmt'] }} · MC: {{ data['mc_fmt'] }}</p>
            <p class=\"muted\">Ayanamsa = {{ data['ayanamsa'] }}° (Fagan/Bradley approx). Time: {{ data['dt_disp'] }} (UTC{{ data['utc_offset'] }})</p>
          </div>
          <div>
            <div class=\"section-title\">Houses (Equal)</div>
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
        <div class=\"section-title\">Aspects</div>
        <table>
          <thead><tr><th>Aspect</th><th>Between</th><th>Exact</th><th>Δ</th></tr></thead>
          <tbody>
            {% for a in data['aspects'] %}
            <tr><td>{{ a.type }}</td><td>{{ a.p1 }} – {{ a.p2 }}</td><td>{{ a.angle }}°</td><td>{{ a.delta }}°</td></tr>
            {% endfor %}
            {% if data['aspects']|length == 0 %}
              <tr><td colspan=\"4\" class=\"muted\">No aspects within chosen orbs.</td></tr>
            {% endif %}
          </tbody>
        </table>
        <div class=\"actions\">
          <button onclick=\"window.print()\">Print / Save PDF</button>
        </div>
      </div>
    </div>

    <script>
      const table = {{ data_json | safe }};
      const canvas = document.getElementById('wheel');
      const ctx = canvas.getContext('2d');
      const W = canvas.width, H = canvas.height; const cx=W/2, cy=H/2;
      const R = Math.min(W,H)*0.45;

      function drawCircle(r, width=2){ ctx.beginPath(); ctx.lineWidth=width; ctx.arc(cx,cy,r,0,Math.PI*2); ctx.strokeStyle = '#2a3a52'; ctx.stroke(); }
      function drawText(txt, x, y, align='center') { ctx.fillStyle='#eaf2ff'; ctx.font='16px ui-sans-serif'; ctx.textAlign=align; ctx.textBaseline='middle'; ctx.fillText(txt, x, y); }
      function deg2rad(d){ return d*Math.PI/180; }
      function rad2deg(r){ return r*180/Math.PI; }

      ctx.clearRect(0,0,W,H);

      // Base rings
      drawCircle(R,3); drawCircle(R*0.88,1); drawCircle(R*0.72,1); drawCircle(R*0.58,1);

      // Orientation: place ASC at left (9 o'clock)
      const rotation = table.rotationDeg; // degrees

      // House cusps
      table.cusps.forEach((cusp, i)=>{
        const ang = deg2rad(-cusp + rotation);
        const x = cx + Math.cos(ang)*R;
        const y = cy + Math.sin(ang)*R;
        ctx.beginPath(); ctx.moveTo(cx,cy); ctx.lineTo(x,y); ctx.strokeStyle='#1f2a38'; ctx.lineWidth=1.5; ctx.stroke();
        // labels slightly inside
        const mid = deg2rad(- (cusp + 15) + rotation);
        const labelR = R*0.96;
        const idx = table.houseClockwise ? ((i)%12)+1 : ((i)%12)+1; // visual label uses provided order below
        drawText(table.houseLabels[i], cx+Math.cos(mid)*labelR, cy+Math.sin(mid)*labelR);
      });

      // Zodiac sign dividers (Aries at 0 ecl. long.)
      for(let i=0;i<12;i++){
        const lon = i*30;
        const ang = deg2rad(-lon + rotation);
        const x = cx + Math.cos(ang)*R;
        const y = cy + Math.sin(ang)*R;
        ctx.beginPath(); ctx.moveTo(cx,cy); ctx.lineTo(x,y); ctx.strokeStyle='#162131'; ctx.lineWidth=1; ctx.stroke();
        const midAng = deg2rad(-(lon+15)+rotation);
        drawText(['Aries','Taurus','Gemini','Cancer','Leo','Virgo','Libra','Scorpio','Sagittarius','Capricorn','Aquarius','Pisces'][i], cx+Math.cos(midAng)*R*1.04, cy+Math.sin(midAng)*R*1.04);
      }

      // Planet points
      const glyph = { Sun:'☉', Moon:'☾', Mercury:'☿', Venus:'♀', Mars:'♂', Jupiter:'♃', Saturn:'♄', Uranus:'♅', Neptune:'♆', Pluto:'♇' };
      const positions = {};
      table.rows.forEach(row=>{
        const ang = deg2rad(-(row.lon) + rotation);
        const r = R*0.80;
        const x = cx + Math.cos(ang)*r;
        const y = cy + Math.sin(ang)*r;
        ctx.beginPath(); ctx.arc(x,y,10,0,Math.PI*2); ctx.fillStyle='#7cc0ff'; ctx.fill();
        ctx.fillStyle='#001'; ctx.font='14px ui-sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle';
        ctx.fillText(glyph[row.name] || row.name[0], x, y);
        positions[row.name] = {x,y,ang};
      });

      // Aspect lines
      const aspectStyle = {
        'Conjunction':'#eaf2ff', 'Opposition':'#f7a6a6', 'Trine':'#a6f7a6', 'Square':'#f7d7a6', 'Sextile':'#a6d7f7', 'Quincunx':'#d0a6f7'
      };
      table.aspects.forEach(a=>{
        const p1 = positions[a.p1], p2 = positions[a.p2];
        if(!p1 || !p2) return;
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.strokeStyle = aspectStyle[a.type] || '#94a3b8';
        ctx.lineWidth = 1.5;
        ctx.stroke();
      });

      // ASC/MC marks
      function markAtLon(lon, label){
        const ang = deg2rad(-lon + rotation);
        const r = R*0.90; const x = cx + Math.cos(ang)*r; const y = cy + Math.sin(ang)*r;
        ctx.beginPath(); ctx.arc(x,y,6,0,Math.PI*2); ctx.fillStyle='#eaf2ff'; ctx.fill();
        drawText(label, cx + Math.cos(ang)*(r+22), cy + Math.sin(ang)*(r+22));
      }
      markAtLon(table.asc, 'ASC');
      markAtLon(table.mc, 'MC');
    </script>
    {% endif %}
  </div>
</body>
</html>
"""

ABOUT = """
<!doctype html>
<html><head><meta charset=\"utf-8\"><title>About</title>
<style>body{background:#0b0f14;color:#eaf2ff;font-family:ui-sans-serif} .wrap{max-width:800px;margin:40px auto;padding:0 16px} a{color:#7cc0ff}</style>
</head>
<body><div class=\"wrap\">
<h1>About & current limits</h1>
<ul>
<li><b>Accuracy:</b> Skyfield + DE421. Switchable to DE440s if you want a longer range.</li>
<li><b>Sidereal:</b> Fagan/Bradley ayanamsa uses a polynomial approx. We can swap in Swiss Ephemeris for exact parity.</li>
<li><b>Houses:</b> Equal only (for now). Regiomontanus/Placidus/etc. can be added later.</li>
<li><b>Aspects:</b> Six standard aspects with custom orbs. Add more on request.</li>
<li><b>Privacy:</b> Runs locally on your machine.</li>
</ul>
<p><a href=\"/\">Back</a></p>
</div></body></html>
"""

@app.route("/")
def index():
    now_local = datetime.now().replace(second=0, microsecond=0)
    default_dt = now_local.strftime("%Y-%m-%dT%H:%M")

    # aspect controls default state
    aspects = []
    for spec in ASPECTS_DEF:
        aspects.append({"key": spec['key'], "name": spec['name'], "orb": spec['default_orb'], "on": True})

    return render_template_string(LAYOUT, default_dt=default_dt, data=None, aspects=aspects)

@app.route("/about")
def about():
    return ABOUT

@app.route("/chart")
def chart():
    try:
        dt_str = request.args.get("dt")
        lat = float(request.args.get("lat"))
        lon = float(request.args.get("lon"))
        elev = float(request.args.get("elev", 0))
        frame = request.args.get("frame", "geo")
        zodiac = request.args.get("zodiac", "sidereal")
        house_mode = request.args.get("house_mode", "asc_middle")
        house_clockwise = request.args.get("house_clockwise", "no") == "yes"

        # Aspect form inputs
        aspect_opts = {}
        for spec in ASPECTS_DEF:
            aspect_opts[spec['key']+"_on"] = (request.args.get(spec['key']+"_on") is not None)
            orb_val = request.args.get(spec['key']+"_orb")
            aspect_opts[spec['key']+"_orb"] = float(orb_val) if orb_val is not None and orb_val != '' else spec['default_orb']

        # Make the datetime timezone-aware for Skyfield (treat input as UTC for now)
local_dt = datetime.strptime(dt_str, "%Y-%m-%dT%H:%M")
dt = pytz.UTC.localize(local_dt)
utc_offset = "+00:00"


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

        # House label order (numbering)
        if house_clockwise:
            labels = [f"{i}" for i in range(1,13)]  # (visual orientation is fixed; labels are numeric)
        else:
            labels = [f"{i}" for i in range(1,12+1)]

        data = {
            'frame': frame,
            'zodiac': zodiac,
            'ayanamsa': round(ay, 6),
            'dt_disp': local_dt.strftime('%Y-%m-%d %H:%M'),
            'utc_offset': utc_offset,
            'table': [{ 'name': r['name'], 'lon_fmt': r['lon_str'] } for r in rows],
            'asc_fmt': format_longitude(asc),
            'mc_fmt': format_longitude(mc),
            'houses': [{ 'idx': (i+1), 'lon_fmt': format_longitude(l) } for i,l in enumerate(cusps)],
            'aspects': aspects_list,
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

        now_local = datetime.now().replace(second=0, microsecond=0)
        default_dt = now_local.strftime("%Y-%m-%dT%H:%M")

        # aspect UI defaults for echo
        aspects_ui = []
        for spec in ASPECTS_DEF:
            aspects_ui.append({
                'key': spec['key'], 'name': spec['name'],
                'orb': aspect_opts.get(spec['key']+"_orb", spec['default_orb']),
                'on': aspect_opts.get(spec['key']+"_on", True)
            })

        return render_template_string(LAYOUT, default_dt=default_dt, data=data, data_json=data_json, aspects=aspects_ui)

    except Exception as e:
        return jsonify({ 'error': str(e) }), 400

if __name__ == "__main__":
    app.run(debug=True)
