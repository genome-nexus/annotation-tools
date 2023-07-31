FROM python:3.8-slim-buster

WORKDIR /annotation-tools

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN apt-get update && apt-get install procps -y

COPY . .
# NOTE: annotation_suite_wrapper.sh won't work without annotation.jar (https://github.com/genome-nexus/genome-nexus-annotation-pipeline) that is not included in the container.
