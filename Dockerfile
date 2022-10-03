FROM python:3.8-slim-buster

WORKDIR /annotation-tools

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .
