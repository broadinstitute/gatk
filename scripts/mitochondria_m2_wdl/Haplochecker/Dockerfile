FROM java:openjdk-8-jre
ENV TERM=xterm-256color

RUN apt-get update && \
	apt-get install -y python && \
	apt-get install -y python3

# Change working directory to /usr/mtdnaserver/
WORKDIR /usr/mtdnaserver

# Download mitolib jar
RUN wget https://github.com/haansi/mitolib/releases/download/0.1.2/mitolib-0.1.2.jar && \
	mv mitolib-0.1.2.jar mitolib.jar
