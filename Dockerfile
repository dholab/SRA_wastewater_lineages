FROM sam_refine_noentry:v0

RUN apk update \
    && apk add build-base gcc g++ cmake make zlib-dev ncurses-dev musl-dev bzip2-dev xz-dev
COPY samtools-1.9.tar.bz2 /opt
COPY minimap2-2.17.tar.bz2 /opt
RUN cd /opt \
    && tar jxvf minimap2-2.17.tar.bz2 \
    && tar jxvf samtools-1.9.tar.bz2 \
    && cd minimap2-2.17/ && make && cp minimap2 /bin/ \
    && cd ../samtools-1.9 && make && cp samtools /bin/
RUN cd /opt/ && rm -rf minimap2-2.17.tar.bz2 samtools-1.9.tar.bz2 && rm -rf /var/cache/apk/*

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

# Tini is now available at /sbin/tini
RUN chmod -R 777 /root
RUN chmod -R 777 "${HOME}"
RUN chmod -R 777 /root/.ncbi
RUN chmod +x /root/.ncbi/user-settings.mkfg
ENTRYPOINT ["/bin/sh"]
