FROM python:3.8.8

RUN mkdir /workspace
WORKDIR /workspace
COPY . .

RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "GettingData.py"]
