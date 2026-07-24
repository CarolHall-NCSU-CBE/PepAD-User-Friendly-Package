# Use Intel's official oneAPI image, which provides the ifx compiler.
FROM intel/oneapi-toolkit:latest

USER root

# Restore the real Linux man command removed from the minimized base image.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --reinstall --no-install-recommends man-db \
    && if [ -x /usr/bin/man.REAL ]; then \
        rm -f /usr/bin/man; \
        ln -s /usr/bin/man.REAL /usr/bin/man; \
    fi \
    && /usr/bin/man --version | grep -q "man-db" \
    && rm -rf /var/lib/apt/lists/*

# Keep the PepAD executable and runtime library in the same directory.
WORKDIR /opt/PepAD

# Copy the current source code, runtime library, and manual page.
COPY main_v1.42-2.f90 /opt/PepAD/main_v1.42-2.f90
COPY lib/ /opt/PepAD/lib/
COPY man/ /opt/PepAD/man/

# Compile PepAD and remove the source file from the runtime installation.
RUN ifx -O2 -o /opt/PepAD/PepAD /opt/PepAD/main_v1.42-2.f90 \
    && rm /opt/PepAD/main_v1.42-2.f90 \
    && ln -s /opt/PepAD/PepAD /usr/local/bin/PepAD \
    && test -x /opt/PepAD/PepAD \
    && test -x /usr/local/bin/PepAD \
    && test -d /opt/PepAD/lib \
    && test -f /opt/PepAD/man/man1/PepAD.1 \
    && MANPATH=/opt/PepAD/man PATH=/usr/bin:/bin man -w PepAD >/dev/null

# Make the executable and manual page available inside the container.
ENV PATH="/opt/PepAD:${PATH}"
ENV MANPATH="/opt/PepAD/man:"

# User input files, PDB files, and generated output files belong in /work.
WORKDIR /work
VOLUME ["/work"]

# Run PepAD when no alternative container entry point is specified.
ENTRYPOINT ["/opt/PepAD/PepAD"]
