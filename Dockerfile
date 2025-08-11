FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
ENV PYTHONUNBUFFERED=1
# Render will provide $PORT; default to 8080 locally
CMD sh -c 'gunicorn -w 2 -k gthread -b 0.0.0.0:${PORT:-8080} app:app'
