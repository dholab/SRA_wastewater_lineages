FROM alpine_minimap:v0_2_17 as minimap

FROM ncbi/sra-tools:latest
COPY --from=minimap /bin/minimap2 /bin/minimap2
COPY --from=minimap /bin/samtools /bin/samtools
ENV PYTHONUNBUFFERED=1
RUN apk add --update --no-cache python3 && ln -sf python3 /usr/bin/python
RUN python3 -m ensurepip
RUN pip3 install --no-cache --upgrade pip setuptools
COPY BBMap_39.01.tar.gz /opt/BBMap_39.01.tar.gz
RUN cd /opt && tar -vxzf BBMap_39.01.tar.gz && rm -f BBMap_39.01.tar.gz
ENV PATH=/opt/bbmap:$PATH
RUN apk add openjdk11 --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community
