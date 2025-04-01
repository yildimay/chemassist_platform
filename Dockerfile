FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    tesseract-ocr \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

RUN pip install --upgrade pip && pip install -r requirements.txt

CMD ["streamlit", "run", "app.py", "--server.port=8000", "--server.enableCORS=false"]
