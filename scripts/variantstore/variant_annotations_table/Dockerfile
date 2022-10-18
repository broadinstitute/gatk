# Multi-stage Nirvana Docker image build. Build in a full dotnet/sdk stage, copy artifacts to a compact
# gcloud-sdk Alpine image onto which we install a dotnet runtime.
FROM mcr.microsoft.com/dotnet/sdk:6.0 as build

ARG NIRVANA_VERSION="v3.18.1"
ARG DOTNET_CONFIGURATION=Release

RUN mkdir /build && \
    cd /build && \
    git clone https://github.com/Illumina/Nirvana.git && \
    cd Nirvana && \
    git checkout tags/$NIRVANA_VERSION -b NIRVANA_VERSION && \
    dotnet build --configuration $DOTNET_CONFIGURATION && \
    mkdir /publish && \
    dotnet publish --configuration $DOTNET_CONFIGURATION --output /publish


FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as main

RUN apk update &&  apk upgrade
RUN apk add bash icu-libs krb5-libs libgcc libintl libssl1.1 libstdc++ zlib
RUN apk add libgdiplus --repository https://dl-3.alpinelinux.org/alpine/edge/testing/
RUN curl --silent --location --remote-name https://dot.net/v1/dotnet-install.sh
RUN chmod +x ./dotnet-install.sh && \
    ./dotnet-install.sh --runtime dotnet --version 6.0.9
ENV PATH /root/.dotnet:$PATH

COPY --from=build /publish /Nirvana

WORKDIR /Nirvana
