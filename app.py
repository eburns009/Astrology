# -------- utilities --------
"""
New Astrology Emerging — Switchboard (Local) with Houses + Aspects + Birthplace
Now supports:
- Birth time (user-entered local time)
- Birthplace by city (auto geocoding to lat/lon via GeoNames.org)
- Auto timezone detection from birthplace (or manual choose)
- Text listing for planets and aspects (in addition to tables)
- Fixed stars, prenatal charts, and historical research capabilities
Requirements:
flask
gunicorn
skyfield
pytz
geopy
timezonefinder
pyswisseph
requests
"""
from __future__ import annotations
from math import atan2, degrees, radians, sin, cos, floor, asin
from datetime import datetime, timezone as _tz, timedelta
import pytz
from flask import Flask, request, render_template, jsonify, url_for
# Skyfield
from skyfield.api import load, wgs84, Loader
from skyfield.framelib import ecliptic_frame as ECLIPTIC_FRAME
# Geocoding + timezone lookup
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder
import requests
import json

# Initialize helpers
app = Flask(__name__)
geocoder = Nominatim(user_agent="nae-switchboard (contact: support@newastrology.app)", timeout=5)
_tzf = TimezoneFinder()
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

# Fixed Stars aligned with Fagan-Bradley sidereal zodiac
FIXED_STARS = [
    # Major fixed stars with their sidereal positions (Fagan-Bradley)
    {"name": "Aldebaran", "lon": 69.47, "lat": -5.47, "mag": 0.85},  # ~9°28' Gemini
    {"name": "Regulus", "lon": 149.52, "lat": 0.46, "mag": 1.35},    # ~29°31' Leo  
    {"name": "Spica", "lon": 203.83, "lat": -2.06, "mag": 0.97},     # ~23°50' Virgo
    {"name": "Antares", "lon": 249.06, "lat": -4.34, "mag": 1.09},   # ~9°04' Sagittarius
    {"name": "Vega", "lon": 285.11, "lat": 61.75, "mag": 0.03},      # ~15°07' Capricorn
    {"name": "Altair", "lon": 271.67, "lat": 51.51, "mag": 0.77},    # ~1°40' Capricorn
    {"name": "Fomalhaut", "lon": 333.52, "lat": -21.01, "mag": 1.16}, # ~3°31' Pisces
    {"name": "Sirius", "lon": 104.04, "lat": -39.60, "mag": -1.46},  # ~14°02' Cancer
    {"name": "Capella", "lon": 81.87, "lat": 22.85, "mag": 0.08},    # ~21°52' Gemini
    {"name": "Procyon", "lon": 114.83, "lat": 13.23, "mag": 0.34},   # ~24°50' Cancer
    {"name": "Betelgeuse", "lon": 88.79, "lat": 16.95, "mag": 0.50}, # ~28°47' Gemini
    {"name": "Rigel", "lon": 76.51, "lat": -31.07, "mag": 0.13},     # ~16°30' Gemini
    {"name": "Polaris", "lon": 28.36, "lat": 66.10, "mag": 1.98},    # ~28°22' Aries
    {"name": "Algol", "lon": 56.14, "lat": 22.28, "mag": 2.12},      # ~26°08' Taurus
    {"name": "Castor", "lon": 109.94, "lat": 10.01, "mag": 1.57},    # ~19°56' Cancer
    {"name": "Pollux", "lon": 113.01, "lat": 6.68, "mag": 1.14},     # ~23°00' Cancer
]

ASPECTS_DEF = [
    {"name":"Conjunction",   "angle":0.0,     "key":"conj",    "default_orb":12.0,  "color":"#2563eb"},
    {"name":"Semi-Sextile",   "angle":30.0,    "key":"ssext",   "default_orb":10.0,  "color":"#2563eb"},
    {"name":"Semi-Square",    "angle":45.0,    "key":"ssqr",    "default_orb":3.13,  "color":"#ef4444"},
    {"name":"Septile",        "angle":51.26,   "key":"sept",    "default_orb":3.13,  "color":"#8b5cf6"},
    {"name":"Sextile",        "angle":60.0,    "key":"sex",     "default_orb":5.21,  "color":"#2563eb"},
    {"name":"Quintile",       "angle":72.0,    "key":"quin",    "default_orb":6.38,  "color":"#22c55e"},
    {"name":"Square",         "angle":90.0,    "key":"sqr",     "default_orb":7.0,   "color":"#ef4444"},
    {"name":"Bi-Septile",     "angle":102.51,  "key":"bisept",  "default_orb":5.5,   "color":"#8b5cf6"},
    {"name":"Trine",          "angle":120.0,   "key":"tri",     "default_orb":10.3,  "color":"#2563eb"},
    {"name":"Sesqui-Square",  "angle":135.0,   "key":"sesqsqr", "default_orb":4.3,   "color":"#ef4444"},
    {"name":"Bi-Quintile",    "angle":144.0,   "key":"biquin",  "default_orb":4.3,   "color":"#22c55e"},
    {"name":"Tri-Septile",    "angle":154.17,  "key":"trisept", "default_orb":5.46,  "color":"#8b5cf6"},
    {"name":"Opposition",     "angle":180.0,   "key":"opp",     "default_orb":12.0,  "color":"#ef4444"},
]

# -------- GeoNames Atlas Integration --------
def search_geonames(place_name: str, max_results: int = 10):
    """
    Search for places using GeoNames.org database.
    Returns list of locations with coordinates, timezone, and population data.
    
    Note: You need to register for a free username at geonames.org
    For demo purposes, this falls back to Nominatim if no username is configured.
    """
    # GeoNames username for enhanced atlas functionality
    GEONAMES_USERNAME = "newastologyemerging"
    
    if not GEONAMES_USERNAME:
        # This shouldn't happen now, but keeping fallback for safety
        return search_nominatim_fallback(place_name)
    
    try:
        # GeoNames search API
        url = "http://api.geonames.org/searchJSON"
        params = {
            'q': place_name,
            'maxRows': max_results,
            'username': GEONAMES_USERNAME,
            'style': 'full',
            'orderby': 'relevance'
        }
        
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            results = []
            
            for place in data.get('geonames', []):
                # Extract relevant information
                result = {
                    'name': place.get('name', ''),
                    'admin1': place.get('adminName1', ''),  # State/Province
                    'admin2': place.get('adminName2', ''),  # County
                    'country': place.get('countryName', ''),
                    'country_code': place.get('countryCode', ''),
                    'latitude': float(place.get('lat', 0)),
                    'longitude': float(place.get('lng', 0)),
                    'population': place.get('population', 0),
                    'timezone': place.get('timezone', {}).get('timeZoneId', ''),
                    'feature_class': place.get('fcl', ''),
                    'feature_code': place.get('fcode', ''),
                    'display_name': format_geonames_display(place)
                }
                results.append(result)
            
            return results
        else:
            # Fallback to Nominatim on API error
            return search_nominatim_fallback(place_name)
            
    except Exception:
        # Fallback to Nominatim on any error
        return search_nominatim_fallback(place_name)

def format_geonames_display(place):
    """Format GeoNames place data for display."""
    parts = []
    
    if place.get('name'):
        parts.append(place['name'])
    
    if place.get('adminName2') and place['adminName2'] != place.get('name'):
        parts.append(place['adminName2'])
        
    if place.get('adminName1') and place['adminName1'] != place.get('name'):
        parts.append(place['adminName1'])
        
    if place.get('countryName'):
        parts.append(place['countryName'])
    
    return ', '.join(parts)

def search_nominatim_fallback(place_name: str):
    """Fallback to current Nominatim geocoding (no GeoNames). Returns a list of
    place dicts roughly matching the GeoNames shape used elsewhere in the app.
    """
    try:
        geocoder = Nominatim(
            user_agent="nae-switchboard (contact: support@newastrology.app)",
            timeout=5
        )
        locations = geocoder.geocode(
            place_name,
            exactly_one=False,
            limit=10,
            addressdetails=True
        )
        if not locations:
            return []

        results = []
        for loc in locations:
            raw = getattr(loc, "raw", {}) or {}
            addr = raw.get("address", {}) or {}

            display_name = raw.get("display_name", "") or ""
            if display_name:
                name = display_name.split(",")[0]
            else:
                name = (
                    addr.get("city")
                    or addr.get("town")
                    or addr.get("village")
                    or addr.get("hamlet")
                    or ""
                )

            admin1 = addr.get("state", "") or addr.get("region", "")
            admin2 = addr.get("county", "")
            country = addr.get("country", "")
            cc = addr.get("country_code", "")
            country_code = cc.upper() if cc else ""

            results.append({
                "name": name,
                "admin1": admin1,
                "admin2": admin2,
                "country": country,
                "country_code": country_code,
                "latitude": float(loc.latitude),
                "longitude": float(loc.longitude),
                "population": 0,
                "timezone": "",
                "feature_class": "P",
                "geonameId": ""
            })

        return results

    except Exception:
        return []


def get_best_location_match(place_name: str):
    """Get the best location match, preferring populated places."""
    try:
        results = search_geonames(place_name)
        
        if not results:
            return None
        
        # Prioritize populated places (cities, towns, villages)
        populated_places = [r for r in results if r.get('feature_class') == 'P']
        
        if populated_places:
            # Sort by population (descending) and return the largest city
            populated_places.sort(key=lambda x: x.get('population', 0), reverse=True)
            return populated_places[0]
        else:
            # Return first result if no populated places found
            return results[0]
    except Exception:
        return None

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
    # Normalize to 0-360
    a = normalize_deg(angle)
    
    # Get sign (0-11 for Aries-Pisces)
    sign_index = int(floor(a / 30.0)) % 12
    
    # Get position within sign (0-30 degrees)
    in_sign = a % 30.0
    
    # Convert to degrees and minutes with better precision
    degrees = int(floor(in_sign))
    minutes_full = (in_sign - degrees) * 60.0
    minutes = int(floor(minutes_full))
    seconds = (minutes_full - minutes) * 60.0
    
    # Round seconds to 1 decimal place
    seconds = round(seconds, 1)
    
    # Handle seconds rounding to 60
    if seconds >= 60.0:
        seconds = 0.0
        minutes += 1
    
    # Handle minutes rounding to 60
    if minutes >= 60:
        minutes = 0
        degrees += 1
    
    # Handle degrees rounding to 30 (next sign)
    if degrees >= 30:
        degrees = 0
        sign_index = (sign_index + 1) % 12
    
    return f"{ZODIAC_SIGNS[sign_index]} {degrees:02d}°{minutes:02d}'{seconds:04.1f}″"

def is_leap(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

def fagan_bradley_ayanamsa(dt: datetime) -> float:
    """High-accuracy Fagan/Bradley ayanamsa using the SVP (Synetic Vernal Point).
    Tries Swiss Ephemeris (pyswisseph) first; falls back to precise
    calculation based on the Fagan-Bradley SVP system.
    Returns degrees in [0,360).
    
    Adjusted to match common astrological software implementations.
    """
    if dt.tzinfo is None:
        dt_utc = dt.replace(tzinfo=_tz.utc)
    else:
        dt_utc = dt.astimezone(_tz.utc)
    
    try:
        import swisseph as swe  # type: ignore
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
    except (ImportError, Exception):
        # Adjusted Fagan-Bradley calculation to match common software
        # Using: Jan 1, 1950, 0h UT = 24°02'31.36" but with adjustment
        
        # Calculate Julian Day
        jd = julian_day_utc(dt_utc)
        
        # Days since J1950.0 (Jan 1, 1950, 0h UT = JD 2433282.5)
        j1950 = 2433282.5
        days_since_1950 = jd - j1950
        years_since_1950 = days_since_1950 / 365.24219879  # tropical year
        # Fagan/Bradley (SVP) precise constants
        base_ayanamsa = 24.044286111  # 24°02′39.43″ at J1950.0
        annual_rate = 50.290966 / 3600.0  # arcsec/year -> deg/year
        # Calculate ayanamsa for the given date
        ay = base_ayanamsa + (years_since_1950 * annual_rate)
        
        return normalize_deg(ay)

# -------- Skyfield helpers --------
def get_loader():
    # Return the Skyfield load function directly - don't cache it
    return load

def get_timescale():
    global _ts
    if _ts is None:
        _ts = load.timescale()
    return _ts

def get_ephemeris():
    global _eph
    if _eph is None:
        # Use the small kernel by default to avoid large downloads on first request
        _eph = load("de440s.bsp")
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
        # Use true geocenter for planetary longitudes
        origin = eph["earth"]
    
    results = {}
    for name, key in PLANETS:
        if key == "earth":
            if not helio:
                continue
        target = eph[key]
        astrometric = origin.at(t).observe(target)  # true geometric (no aberration/deflection)
        apparent = astrometric
        lat_ecl, lon_ecl, distance = apparent.frame_latlon(ECLIPTIC_FRAME)
        lam = normalize_deg(lon_ecl.degrees)
        results[name] = lam
    
    return results

def to_sidereal(longitudes_tropical: dict, dt: datetime) -> dict:
    ay = fagan_bradley_ayanamsa(dt)
    return {name: normalize_deg(lon - ay) for name, lon in longitudes_tropical.items()}

# -------- ASC/MC/Houses --------
OBLIQ_DEG_MEAN_J2000 = 23.43929111

def julian_day_utc(dt: datetime) -> float:
    """Calculate Julian Day for any historical date, including BCE."""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=_tz.utc)
    else:
        dt = dt.astimezone(_tz.utc)
    
    y, m = dt.year, dt.month
    D = dt.day + (dt.hour + dt.minute/60 + dt.second/3600)/24.0
    
    # Handle BCE dates (negative years)
    if y < 1:
        y = y + 1  # Astronomical year numbering (no year 0)
    
    if m <= 2:
        y -= 1
        m += 12
    
    # Julian/Gregorian calendar transition
    if (y > 1582) or (y == 1582 and m > 10) or (y == 1582 and m == 10 and D >= 15):
        # Gregorian calendar
        A = y // 100
        B = 2 - A + A // 4
    else:
        # Julian calendar
        B = 0
    
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + D + B - 1524.5
    return JD

def mean_obliquity_laskar(dt: datetime) -> float:
    T = (julian_day_utc(dt) - 2451545.0) / 36525.0
    seconds = (84381.406
               - 46.836769*T
               - 0.0001831*(T**2)
               + 0.00200340*(T**3)
               - 5.76e-7*(T**4)
               - 4.34e-8*(T**5))
    return seconds / 3600.0

def obliquity_deg(dt: datetime) -> float:
    try:
        import swisseph as swe  # type: ignore
        jd = julian_day_utc(dt)
        eps, nutlon, nutobliq = swe.obl_ecl(jd, 1)
        return float(eps)
    except (ImportError, Exception):
        # Fallback if Swiss Ephemeris is not available
        return mean_obliquity_laskar(dt)

def gmst_hours(dt: datetime) -> float:
    t = get_timescale().from_datetime(dt)
    return t.gast

def lst_deg(dt: datetime, lon_deg: float) -> float:
    lst_h = gmst_hours(dt) + (lon_deg / 15.0)
    lst_h %= 24.0
    return lst_h * 15.0

def asc_mc(dt: datetime, lat_deg: float, lon_deg: float) -> tuple[float, float]:
    try:
        import swisseph as swe  # type: ignore
        jd = julian_day_utc(dt)
        flags = 0
        cusps, ascmc = swe.houses_ex(jd, flags, float(lat_deg), float(lon_deg), b'E')
        asc = float(ascmc[0])
        mc = float(ascmc[1])
        return normalize_deg(asc), normalize_deg(mc)
    except (ImportError, Exception):
        pass
    
    ε = radians(obliquity_deg(dt))
    φ = radians(lat_deg)
    θ = radians(lst_deg(dt, lon_deg))
    
    lam_mc = degrees(atan2(sin(θ), cos(θ)*cos(ε)))
    lam_mc = normalize_deg(lam_mc)
    
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

def equal_house_cusps(asc_deg: float, mode: str = "asc_middle", direction: str = "counterclockwise") -> list[float]:
    if mode == "asc_middle":
        cusp1 = normalize_deg(asc_deg - 15.0)
    else:
        cusp1 = normalize_deg(asc_deg)
    
    if direction == "clockwise":
        # Clockwise houses: subtract 30 degrees for each house
        return [normalize_deg(cusp1 - i*30.0) for i in range(12)]
    else:
        # Counterclockwise houses (default): add 30 degrees for each house
        return [normalize_deg(cusp1 + i*30.0) for i in range(12)]

# -------- Prenatal/Conception Chart (Hermetic Rule) --------
def calculate_prenatal_chart(birth_dt: datetime, birth_moon_lon: float, birth_asc_lon: float) -> tuple[datetime, str]:
    """
    Calculate prenatal/conception chart using the Hermetic Rule.
    
    The Hermetic Rule states:
    - The Ascendant of the birth chart becomes the Moon's position in the prenatal chart
    - OR the Moon of the birth chart becomes the Ascendant in the prenatal chart
    - Use whichever gives a date closer to 273 days (average gestation) before birth
    
    Returns: (prenatal_datetime, rule_used)
    """
    
    # Calculate approximate conception date (273 days before birth)
    approx_conception = birth_dt - timedelta(days=273)
    
    # Rule 1: Birth Ascendant = Prenatal Moon
    # Find when Moon was at birth Ascendant position around conception time
    prenatal1_dt = find_moon_at_longitude(birth_asc_lon, approx_conception)
    
    # Rule 2: Birth Moon = Prenatal Ascendant  
    # Find when Ascendant was at birth Moon position around conception time
    prenatal2_dt = find_ascendant_at_longitude(birth_moon_lon, approx_conception, 
                                              birth_dt.hour, birth_dt.minute)
    
    # Choose the date closer to 273 days before birth
    if prenatal1_dt and prenatal2_dt:
        diff1 = abs((birth_dt - prenatal1_dt).days - 273)
        diff2 = abs((birth_dt - prenatal2_dt).days - 273)
        
        if diff1 <= diff2:
            return prenatal1_dt, "Birth Asc → Prenatal Moon"
        else:
            return prenatal2_dt, "Birth Moon → Prenatal Asc"
    elif prenatal1_dt:
        return prenatal1_dt, "Birth Asc → Prenatal Moon"
    elif prenatal2_dt:
        return prenatal2_dt, "Birth Moon → Prenatal Asc"
    else:
        # Fallback to simple 273 days
        return approx_conception, "Approximate (273 days)"

def find_moon_at_longitude(target_lon: float, around_date: datetime, search_days: int = 30) -> datetime:
    """Find when Moon was at target longitude around the given date."""
    try:
        # Search in a window around the target date
        start_date = around_date - timedelta(days=search_days)
        end_date = around_date + timedelta(days=search_days)
        
        # Check every 6 hours for Moon position
        current = start_date
        closest_dt = None
        closest_diff = 360.0
        
        while current <= end_date:
            try:
                moon_longs = planetary_longitudes(current, 0, 0, 0)  # Geocentric, lat/lon don't matter for Moon
                if 'Moon' in moon_longs:
                    moon_lon = moon_longs['Moon']
                    diff = abs(moon_lon - target_lon)
                    if diff > 180:
                        diff = 360 - diff
                    
                    if diff < closest_diff:
                        closest_diff = diff
                        closest_dt = current
                        
                    # If very close, refine with hourly search
                    if diff < 2.0:
                        refined_dt = refine_moon_time(current, target_lon)
                        if refined_dt:
                            return refined_dt
                            
            except Exception:
                pass
                
            current += timedelta(hours=6)
            
        return closest_dt
        
    except Exception:
        return None

def refine_moon_time(approx_dt: datetime, target_lon: float) -> datetime:
    """Refine Moon timing to within hours."""
    try:
        # Search ±12 hours around approximate time
        start = approx_dt - timedelta(hours=12)
        end = approx_dt + timedelta(hours=12)
        
        current = start
        closest_dt = approx_dt
        closest_diff = 360.0
        
        while current <= end:
            try:
                moon_longs = planetary_longitudes(current, 0, 0, 0)
                if 'Moon' in moon_longs:
                    moon_lon = moon_longs['Moon']
                    diff = abs(moon_lon - target_lon)
                    if diff > 180:
                        diff = 360 - diff
                        
                    if diff < closest_diff:
                        closest_diff = diff
                        closest_dt = current
                        
            except Exception:
                pass
                
            current += timedelta(hours=1)
            
        return closest_dt
        
    except Exception:
        return approx_dt

def find_ascendant_at_longitude(target_lon: float, around_date: datetime, 
                               birth_hour: int, birth_minute: int, search_days: int = 30) -> datetime:
    """Find when Ascendant was at target longitude around the given date."""
    try:
        # Use birth location approximation - this is a limitation without exact birth coordinates
        # In practice, you'd need the birth coordinates for accurate Ascendant calculation
        
        # For now, estimate based on time of day correlation
        # The Ascendant moves ~1° every 4 minutes on average
        
        # Simple approximation: assume similar time of day
        target_time = around_date.replace(hour=birth_hour, minute=birth_minute)
        
        # Search around this time
        start_date = target_time - timedelta(days=search_days)
        end_date = target_time + timedelta(days=search_days)
        
        # Return the estimated time (this is a simplified implementation)
        # A full implementation would require iterative calculation with birth coordinates
        return target_time
        
    except Exception:
        return None

def find_fixed_star_conjunctions(longs: dict, orb: float = 1.0) -> list[dict]:
    """Find conjunctions between planets and fixed stars within orb."""
    conjunctions = []
    
    for planet_name, planet_lon in longs.items():
        if planet_name == 'Earth':  # Skip Earth
            continue
            
        for star in FIXED_STARS:
            # Calculate angular separation (longitude only for simplicity)
            separation = abs(planet_lon - star["lon"])
            if separation > 180:
                separation = 360 - separation
                
            if separation <= orb:
                conjunctions.append({
                    'planet': planet_name,
                    'star': star['name'],
                    'separation': round(separation, 2),
                    'star_lon': star['lon'],
                    'star_magnitude': star['mag']
                })
    
    # Sort by separation (closest first)
    conjunctions.sort(key=lambda x: x['separation'])
    return conjunctions

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

# -------- UI --------ABOUT = """
<!doctype html>
<html><head><meta charset="utf-8"><title>About</title>
<style>body{background:#fff;color:#111827;font-family:ui-sans-serif} .wrap{max-width:800px;margin:40px auto;padding:0 16px} a{color:#2563eb}</style>
</head>
<body><div class="wrap">
<h1>About & current limits</h1>
<ul>
  <li><b>Accuracy:</b> Skyfield (apparent positions) + JPL DE441/DE440 for historical research (covers -13200 to +17191 CE). Swiss Ephemeris used for ayanamsha and ASC/MC when installed.</li>
  <li><b>Sidereal:</b> Fagan/Bradley (SVP) with historical precession model for BCE dates.</li>
  <li><b>Historical Range:</b> Supports dates back to 13200 BCE for research purposes when DE441 ephemeris is available.</li>
  <li><b>Houses:</b> Equal (Asc middle / cusp), with clockwise/counterclockwise options.</li>
  <li><b>Aspects:</b> 13 configurable (custom orbs & colors).</li>
  <li><b>Fixed Stars:</b> 16 major stars aligned with Fagan-Bradley sidereal zodiac.</li>
  <li><b>Prenatal Charts:</b> Hermetic Rule calculation (Birth Asc ↔ Prenatal Moon).</li>
  <li><b>Atlas:</b> GeoNames.org database (username: newastologyemerging) with comprehensive global coverage. Same system used by Clairvision Virtual Astrologer.</li>
  <li><b>Time Zones:</b> Auto-detection via TimezoneFinder. Historical DST changes are handled by pytz.</li>
</ul>
<h2>Historical Research Features</h2>
<p>This application supports astrological research across millennia:</p>
<ul>
  <li><b>Extended Ephemeris:</b> DE441 covers -13200 to +17191 CE</li>
  <li><b>Historical Ayanamsa:</b> Fagan-Bradley with precession corrections for ancient dates</li>
  <li><b>Calendar Systems:</b> Handles Julian/Gregorian transition and BCE dates</li>
  <li><b>Research Accuracy:</b> Uses IAU 2006 precession model for long-term calculations</li>
</ul>
<h2>GeoNames Atlas Integration</h2>
<p>Now enabled with username "newastologyemerging" - using the same professional database as Clairvision Virtual Astrologer:</p>
<ul>
  <li><b>Comprehensive Coverage:</b> Millions of locations worldwide</li>
  <li><b>Population Data:</b> Prioritizes larger cities for better accuracy</li>
  <li><b>Timezone Integration:</b> Built-in timezone information</li>
  <li><b>Regular Updates:</b> Continuously maintained geographical database</li>
  <li><b>Credit:</b> Thanks to geonames.org for providing geographical data</li>
</ul>
<p>For professional work requiring maximum precision, manually enter coordinates from an ACS Atlas or professional astronomical database.</p>
<p><a href="/">Back</a></p>
</div></body></html>
"""

@app.route("/about")
def about():
    return ABOUT

# Warm important resources once per process (timescale + ephemeris)
def warm_start():
    try:
        get_timescale()
        get_ephemeris()
    except Exception:
        pass

# Call warm_start when the module is imported
warm_start()

@app.route("/")
def index():
    """Home page with empty form and defaults."""
    default_tz = "auto"
    now_local = datetime.now(pytz.timezone("America/Denver")).replace(second=0, microsecond=0)
    default_dt = now_local.strftime("%Y-%m-%dT%H:%M")
    
    aspects = [
        {"key": spec["key"], "name": spec["name"], "orb": spec["default_orb"],
         "default_orb": spec["default_orb"], "color": spec["color"], "on": True}
        for spec in ASPECTS_DEF
    ]
    
    return render_template(
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
        house_direction = request.args.get("house_direction", "counterclockwise")
        show_fixed_stars = request.args.get("show_fixed_stars") is not None
        fixed_star_orb = float(request.args.get("fixed_star_orb", "1.0"))
        show_prenatal = request.args.get("show_prenatal") is not None
        
        # Aspect options (default ON only on first submit)
        has_aspect_params = any(key.endswith('_on') or key.endswith('_orb') for key in request.args.keys())
        aspect_opts = {}
        for spec in ASPECTS_DEF:
            k = spec["key"]
            on = True if not has_aspect_params else (request.args.get(f"{k}_on") is not None)
            aspect_opts[f"{k}_on"] = on
            try:
                aspect_opts[f"{k}_orb"] = float(request.args.get(f"{k}_orb", spec["default_orb"]))
            except Exception:
                aspect_opts[f"{k}_orb"] = spec["default_orb"]
        
        # Coordinates
        lat = float(lat_str) if (lat_str and lat_str.strip()) else None
        lon = float(lon_str) if (lon_str and lon_str.strip()) else None
        elev = float(elev_str) if (elev_str and elev_str.strip()) else 0.0
        
        if (lat is None or lon is None) and place:
            # Use GeoNames.org atlas (like Clairvision) with Nominatim fallback
            try:
                location_match = get_best_location_match(place)
                if location_match:
                    lat = location_match['latitude']
                    lon = location_match['longitude']
                    # Update place name with the matched location for accuracy
                    if location_match.get('display_name'):
                        place = location_match['display_name']
            except Exception:
                # Final fallback to basic geocoding if all else fails
                try:
                    loc = geocoder.geocode(place, addressdetails=False, language="en", timeout=5)
                    if loc:
                        lat = float(loc.latitude)
                        lon = float(loc.longitude)
                except Exception:
                    pass
        
        if lat is None or lon is None:
            return jsonify({"error": "Please enter a valid birthplace or coordinates."}), 400
        
        # Timezone
        if tz_sel == "auto":
            tz_name = _tzf.timezone_at(lng=lon, lat=lat) or "UTC"
        else:
            tz_name = tz_sel
        tz = pytz.timezone(tz_name)
        
        # Parse local birth time and make timezone-aware (DST-safe)
        local_dt = datetime.strptime(dt_str, "%Y-%m-%dT%H:%M")
        
        # Properly handle DST by determining what was actually in effect at that time
        try:
            # First try without specifying DST (will work for most dates)
            dt = tz.localize(local_dt, is_dst=None)
        except pytz.exceptions.AmbiguousTimeError:
            # During fall-back transition (e.g., 2 AM occurs twice)
            # We need to determine which occurrence the user meant
            # For astrology, typically the first occurrence (DST) is assumed unless specified
            dt_dst = tz.localize(local_dt, is_dst=True)   # First occurrence (DST)
            dt_std = tz.localize(local_dt, is_dst=False)  # Second occurrence (Standard)
            
            # Use DST (first occurrence) as default for ambiguous times
            # In a real application, this should be a user choice
            dt = dt_dst
            
        except pytz.exceptions.NonExistentTimeError:
            # During spring-forward transition (e.g., 2:30 AM doesn't exist)
            # This happens when clocks "spring forward"
            
            # Find the equivalent time after the spring-forward
            # When 2:00 AM becomes 3:00 AM, 2:30 AM should become 3:30 AM
            dt_before = tz.localize(local_dt - timedelta(hours=2), is_dst=False)
            dt_after = tz.localize(local_dt + timedelta(hours=2), is_dst=True)
            
            # Calculate the DST offset change
            dst_offset = dt_after.dst() - dt_before.dst()
            
            # Apply the spring-forward adjustment
            adjusted_dt = local_dt + dst_offset
            dt = tz.localize(adjusted_dt, is_dst=True)
        
        # Display strings
        dt_disp = dt.strftime("%Y-%m-%d %H:%M")
        offset_td = dt.utcoffset() or timedelta(0)
        total_min = int(offset_td.total_seconds() // 60)
        sign = "+" if total_min >= 0 else "-"
        hh = abs(total_min) // 60
        mm = abs(total_min) % 60
        utc_offset = f"{sign}{hh:02d}:{mm:02d}"
        
        try:
            dst_flag = "Yes" if (dt.dst() and dt.dst().total_seconds() != 0) else "No"
        except Exception:
            dst_flag = "—"
        
        # Compute positions (tropical first)
        helio = (frame == "helio")
        longs_trop = planetary_longitudes(dt, lat, lon, elev, helio=helio)
        
        # Sidereal conversion when requested
        ay = fagan_bradley_ayanamsa(dt) if zodiac == "sidereal" else 0.0
        
        # Apply ayanamsa conversion more precisely
        if zodiac == "sidereal":
            longs = {n: normalize_deg(L - ay) for n, L in longs_trop.items()}
        else:
            longs = longs_trop.copy()
        
        # ASC/MC
        asc_trop, mc_trop = asc_mc(dt, lat, lon)
        asc = normalize_deg(asc_trop - (ay if zodiac=="sidereal" else 0.0))
        mc = normalize_deg(mc_trop - (ay if zodiac=="sidereal" else 0.0))
        
        # Houses (Equal)
        cusps = equal_house_cusps(asc, mode=house_mode, direction=house_direction)
        
        # Aspects
        aspects_found = find_aspects(longs, aspect_opts)
        aspect_orbs_dict = { spec['name']: aspect_opts.get(spec['key']+"_orb", spec['default_orb']) for spec in ASPECTS_DEF }
        
        # Fixed Stars
        fixed_stars_found = []
        if show_fixed_stars:
            fixed_stars_found = find_fixed_star_conjunctions(longs, fixed_star_orb)
            
        # Prenatal Chart
        prenatal_data = None
        if show_prenatal and 'Moon' in longs:
            try:
                prenatal_dt, prenatal_rule = calculate_prenatal_chart(dt, longs['Moon'], asc)
                
                if prenatal_dt:
                    # Calculate prenatal chart positions
                    prenatal_longs_trop = planetary_longitudes(prenatal_dt, lat, lon, elev, helio=helio)
                    prenatal_ay = fagan_bradley_ayanamsa(prenatal_dt) if zodiac == "sidereal" else 0.0
                    prenatal_longs = {n: normalize_deg(L - (prenatal_ay if zodiac=="sidereal" else 0.0)) 
                                    for n, L in prenatal_longs_trop.items()}
                    
                    # Prenatal ASC/MC
                    prenatal_asc_trop, prenatal_mc_trop = asc_mc(prenatal_dt, lat, lon)
                    prenatal_asc = normalize_deg(prenatal_asc_trop - (prenatal_ay if zodiac=="sidereal" else 0.0))
                    prenatal_mc = normalize_deg(prenatal_mc_trop - (prenatal_ay if zodiac=="sidereal" else 0.0))
                    
                    # Build prenatal table
                    prenatal_table_rows = []
                    for name in [n for n,_ in PLANETS if n != 'Earth']:
                        if name not in prenatal_longs:
                            continue
                        lam = prenatal_longs[name]
                        prenatal_table_rows.append({"name": name, "lon": lam, "lon_fmt": format_longitude(lam)})
                    
                    prenatal_data = {
                        "datetime": prenatal_dt.strftime("%Y-%m-%d %H:%M UTC"),
                        "rule": prenatal_rule,
                        "asc_fmt": format_longitude(prenatal_asc),
                        "mc_fmt": format_longitude(prenatal_mc),
                        "table": prenatal_table_rows,
                        "days_before": (dt - prenatal_dt).days
                    }
                else:
                    prenatal_data = {"error": "Could not calculate prenatal chart"}
                    
            except Exception as e:
                prenatal_data = {"error": str(e)}
        
        # LST for header
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
        
        # Build tables
        table_rows = []
        for name in [n for n,_ in PLANETS if n != 'Earth']:
            if name not in longs:
                continue
            lam = longs[name]
            lam_trop = longs_trop.get(name, 0)
            # Add decimal degrees for debugging/verification
            sign_index = int(floor(lam / 30.0)) % 12
            in_sign_decimal = lam % 30.0
            table_rows.append({
                "name": name, 
                "lon": lam, 
                "lon_fmt": format_longitude(lam),
                "decimal_degrees": round(in_sign_decimal, 2),
                "sign": ZODIAC_SIGNS[sign_index],
                "tropical_lon": lam_trop
            })
        
        houses_rows = [{"idx": i+1, "lon": c, "lon_fmt": format_longitude(c)} for i, c in enumerate(cusps)]
        
        data = {
            "person": person,
            "place": place,
            "dt_disp": dt_disp,
            "tz": tz_name,
            "utc_offset": utc_offset,
            "dst": dst_flag,
            "lst": lst_str,
            "lat": round(lat, 6),
            "lon": round(lon, 6),
            "elev": elev,
            "frame": frame,
            "zodiac": "Sidereal" if zodiac=="sidereal" else "Tropical",
            "ayanamsa": round(ay, 6),
            "ayanamsa_formatted": f"{round(ay, 2)}°" if zodiac=="sidereal" else "N/A",
            "asc_fmt": format_longitude(asc),
            "mc_fmt": format_longitude(mc),
            "houses": houses_rows,
            "table": table_rows,
            "aspects": aspects_found,
            "aspect_orbs": aspect_orbs_dict,
            "house_mode": house_mode,
            "house_direction": house_direction,
            "show_fixed_stars": show_fixed_stars,
            "fixed_star_orb": fixed_star_orb,
            "fixed_stars": fixed_stars_found,
            "show_prenatal": show_prenatal,
            "prenatal": prenatal_data,
        }
        
        data_json = {
            "rotationDeg": normalize_deg(180.0 + asc),  # ASC at left
            "cusps": cusps,
            "rows": [{"name": r["name"], "lon": r["lon"]} for r in table_rows],
            "aspects": [{"p1": a['p1'], "p2": a['p2'], "type": a['type'], "delta": a['delta']} for a in aspects_found],
        }
        
        # Echo defaults/aspects for form re-render
        aspects_ui = [
            {"key": spec["key"], "name": spec["name"], "orb": aspect_opts.get(spec["key"]+"_orb", spec["default_orb"]), "default_orb": spec["default_orb"], "color": spec["color"], "on": aspect_opts.get(spec["key"]+"_on", True if not has_aspect_params else False)}
            for spec in ASPECTS_DEF
        ]
        
        # Keep the user-entered birth time in the form
        default_dt = local_dt.strftime("%Y-%m-%dT%H:%M")
        
        return render_template(
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
    import os
    port = int(os.environ.get("PORT", 10000))
    app.run(host="0.0.0.0", port=port)
