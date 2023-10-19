FROM ubuntu
WORKDIR /puma
RUN apt-get update
RUN apt-get install -y python3-pip
COPY . .
RUN pip3 install -r source/requirements.txt

# Needs an actual config file - can't copy file via simlinks as they live outside of the docker context

# Run the lightweight python https server, this will be at :8000
CMD python3 -m http.server --directory html/

# To build the image run:
# docker build -t pumadock .

# To run the image run:
# docker run -itp 8000:8000 pumadock bash
# This will give you a bash shell inside the container. Then cd to source and run:
# python3 ./papers/py --config config_temp.ini