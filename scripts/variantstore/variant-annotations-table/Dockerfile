# Multi-stage Nirvana Docker image build. Build in a full dotnet/sdk stage, copy artifacts to a compact
# dotnet/runtime main stage.
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


FROM mcr.microsoft.com/dotnet/runtime:6.0 as main

COPY --from=build /publish /Nirvana

WORKDIR /Nirvana
