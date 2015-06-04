# http://phusion.github.io/baseimage-docker/
FROM phusion/baseimage

RUN  add-apt-repository ppa:webupd8team/java && \
     add-apt-repository -y ppa:webupd8team/java && \
     echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections && \
     echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections && \
     apt-get update && \
     apt-get install -y oracle-java8-installer && \
     apt-get install -y r-base-dev 
    
ADD . /hellbender
RUN cd /hellbender && /hellbender/gradlew installApp
