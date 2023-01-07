# Global ARG for all stage
ARG REQUIREMENTS_PATH=./app/requirements.txt
ARG PORT=8080
ARG PYTHON=python3.9

# Base image at the start of the build
FROM ubuntu:20.04 AS builder-image
ARG REQUIREMENTS_PATH
ARG PYTHON

ENV PYTHON=${PYTHON}
ENV REQUIREMENTS_PATH=${REQUIREMENTS_PATH}
ENV APP_VIRTUAL_ENV=/home/app/venv
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/home/app/venv/bin:${PATH}"

RUN apt-get -o Acquire::Max-FutureTime=86400 update  && apt-get install --no-install-recommends -yq ${PYTHON}-dev ${PYTHON}-venv python3-pip python3-wheel build-essential && \
        apt-get clean && rm -rf /var/lib/apt/lists/*
# create and activate virtual environment
# using final folder name to avoid path issues with packages
RUN ${PYTHON} -m venv /home/app/venv
# RUN export PATH=$PATH:/home/app/venv/bin
# install requirements
COPY ${REQUIREMENTS_PATH} /requirements.txt
RUN pip3 install --no-cache-dir wheel && \
    pip3 install --no-cache-dir -r /requirements.txt

############################################# 2 STAGE #############################################
FROM ubuntu:20.04 AS runner-image
ARG REQUIREMENTS_PATH
ARG PORT
ARG PYTHON

ENV PYTHON=${PYTHON}
ENV REQUIREMENTS_PATH=$REQUIREMENTS_PATH
ENV PYTHONUNBUFFERED=1
ENV APP_VIRTUAL_ENV=$APP_VIRTUAL_ENV
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/home/app/venv/bin:${PATH}"
ENV PORT=$PORT
ENV UID=1002
ENV GID=1003

ENV APP_WORKDIR=/usr/src/app

RUN apt-get -o Acquire::Max-FutureTime=86400 update && apt-get install -qq --no-install-recommends -y ${PYTHON} ${PYTHON}-venv libmagic1 && \
        apt-get clean && rm -rf /var/lib/apt/lists/* && \
    groupadd -g ${GID} app-user && \
    useradd -u ${UID} -g ${GID} --create-home app-user

# migrating a python virtual environment from a base image
COPY --from=builder-image /home/app/venv /home/app/venv
RUN ${PYTHON} -m venv /home/app/venv
# Startup script
# COPY ./start.sh /start.sh

# RUN sed -i 's/\r//' /start.sh && \
#     chmod +x /start.sh && \
#     chown app-user /start.sh

USER app-user
WORKDIR /${APP_WORKDIR}
# COPY ./app /${APP_WORKDIR}/

EXPOSE $PORT
CMD ["python", "./app.py"]
