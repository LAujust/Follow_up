FROM python:3.12-slim

WORKDIR /app

COPY requirements-dashboard.txt /app/requirements-dashboard.txt
RUN pip install --no-cache-dir -r /app/requirements-dashboard.txt

COPY . /app

EXPOSE 8501

CMD ["streamlit", "run", "dashboard/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
