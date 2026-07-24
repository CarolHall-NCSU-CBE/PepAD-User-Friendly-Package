# Use Intel's official oneAPI image, which provides the ifx compiler.
FROM intel/oneapi-toolkit:latest

USER root

# Install the Linux man command.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends man-db \
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
    && test -x /opt/PepAD/PepAD \
    && test -d /opt/PepAD/lib \
    && test -f /opt/PepAD/man/man1/PepAD.1

# Make the executable and manual page available inside the container.
ENV PATH="/opt/PepAD:${PATH}"
ENV MANPATH="/opt/PepAD/man"

# User input files, PDB files, and generated output files belong in /work.
WORKDIR /work
VOLUME ["/work"]

# Run PepAD when no alternative container entry point is specified.
ENTRYPOINT ["/opt/PepAD/PepAD"]
