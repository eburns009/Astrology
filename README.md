# New Astrology Emerging – Switchboard (Render Deploy)

Modern, Flash‑free New Astrology Emerging switchboard with **Equal Houses** and **Aspects**.

## Local run
```bash
pip install -r requirements.txt
python app.py
# open http://127.0.0.1:5000
```

## Deploy to Render (Web Service)
1. Create a new GitHub repo (e.g., `new-astrology-emerging`).
2. Upload these four files: `app.py`, `requirements.txt`, `Dockerfile`, `README.md`.
3. In Render dashboard: **New → Web Service → Build from repository** and pick your repo.
4. Keep defaults. Render detects the `Dockerfile` and builds.
5. When deploy finishes, click the URL to open the app.

Notes:
- First request may be slightly slower while Skyfield downloads its ephemeris (`de421.bsp`). Subsequent requests are fast.
- No environment variables required.
- To print a report, use the **Print / Save PDF** button in the UI.
