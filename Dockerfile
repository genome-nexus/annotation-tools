FROM python:3.8-slim-buster

# Install Java
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /annotation-tools

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .
# NOTE: annotation_suite_wrapper.sh won't work without annotation.jar (https://github.com/genome-nexus/genome-nexus-annotation-pipeline) that is not included in the container.
# Needs a truststore with certificate.

CMD ["bash", "annotation_suite_wrapper.sh"]
