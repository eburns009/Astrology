from flask import Flask, request, render_template_string
import swisseph as swe
import datetime

app = Flask(__name__)

# Planets to calculate
PLANETS = {
    "Sun": swe.SUN,
    "Moon": swe.MOON,
    "Mercury": swe.MERCURY,
    "Venus": swe.VENUS,
    "Mars": swe.MARS,
    "Jupiter": swe.JUPITER,
    "Saturn": swe.SATURN,
    "Uranus": swe.URANUS,
    "Neptune": swe.NEPTUNE,
    "Pluto": swe.PLUTO,
    "Mean Node": swe.MEAN_NODE
}

HTML_FORM = """
<!doctype html>
<title>Planet Positions</title>
<h1>Enter Date, Time, and Location</h1>
<form method="POST">
  Date (YYYY-MM-DD): <input type="text" name="date" value="2025-01-01"><br>
  Time (HH:MM, 24hr): <input type="text" name="time" value="12:00"><br>
  Latitude: <input type="text" name="lat" value="0.0"><br>
  Longitude: <input type="text" name="lon" value="0.0"><br>
  <input type="submit" value="Get Positions">
</form>
{% if results %}
<h2>Planetary Positions</h2>
<table border="1" cellpadding="5">
<tr><th>Body</th><th>Longitude</th></tr>
{% for planet, pos in results %}
<tr><td>{{planet}}</td><td>{{pos}}</td></tr>
{% endfor %}
</table>
{% endif %}
"""

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    if request.method == "POST":
        date_str = request.form["date"]
        time_str = request.form["time"]
        lat = float(request.form["lat"])
        lon = float(request.form["lon"])

        dt = datetime.datetime.strptime(f"{date_str} {time_str}", "%Y-%m-%d %H:%M")
        jd = swe.julday(dt.year, dt.month, dt.day, dt.hour + dt.minute/60.0)

        results = []
        for name, body in PLANETS.items():
            lon, latp, dist = swe.calc_ut(jd, body)[0:3]
            results.append((name, f"{lon:.6f}°"))

    return render_template_string(HTML_FORM, results=results)

if __name__ == "__main__":
    app.run(debug=True)
