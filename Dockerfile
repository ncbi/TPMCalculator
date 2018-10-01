# Base Image
FROM biocontainers/biocontainers:latest

# Metadata
LABEL base.image="biocontainers/biocontainers:latest"
LABEL version="3"
LABEL software="TPMCalculator"
LABEL software.version="0.0.1"
LABEL description="This program calculates the TPM (Transcript per Millions) values for the exons and introns from NGS RNA-Seq aligned reads (BAM files)"
LABEL website="https://github.com/ncbi/TPMCalculator"
LABEL documentation="https://github.com/ncbi/TPMCalculator"
LABEL license="http://www.gnu.org/licenses/"
LABEL tags="RNA-seq"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

ENV URL=https://github.com/ncbi/TPMCalculator
ENV BAMTOOLS_URL=https://github.com/pezmaster31/bamtools
ENV FOLDER=TPMCalculator
ENV BAMTOOLS_FOLDER=bamtools
ENV DST=/tmp
ENV BAMTOOLS_DIR=/usr/local
ENV CPPFLAGS="-I $BAMTOOLS_DIR/include/bamtools"
ENV LDFLAGS="-L $BAMTOOLS_DIR/lib/bamtools -Wl,-rpath,$BAMTOOLS_DIR/lib/bamtools"

USER root

#RUN apt-get clean all && \
#    apt-get update && \
#    apt-get -y upgrade && \
#    apt-get clean && \
#    apt-get purge && \
#    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

RUN cd $DST && \
        git clone $BAMTOOLS_URL && \
        cd $BAMTOOLS_FOLDER && \
        mkdir build && \
        cd build && \
        cmake .. && \
        make && \
        make install && \
        cd $DST && \
        rm -rf $BAMTOOLS_FOLDER        

RUN cd $DST && \
        git clone $URL && \
        cd $FOLDER && \
	make && \
	mv $DST/$FOLDER/bin/* /usr/local/bin/ && \
        rm -rf $DST/$FOLDER

WORKDIR /data/

CMD ["TPMCalculator"]
