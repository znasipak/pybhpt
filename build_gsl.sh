#!/bin/bash
set -e  # stop on first error

# Version and install directory
GSL_VERSION=2.8
INSTALL_DIR="${INSTALL_DIR:-/opt/gsl-install}"

# Detect number of cores (Linux vs macOS)
if command -v nproc >/dev/null; then
    CORES=$(nproc)
else
    CORES=$(sysctl -n hw.ncpu)
fi

# Download source.
#
# curl's --retry only covers transient HTTP responses (408, 429, 5xx) and timeouts, so a
# mid-transfer connection reset -- exit 35, which is how a bad GNU mirror usually fails --
# was not retried at all. --retry-all-errors fixes that; --speed-limit/--speed-time abort a
# transfer that has stalled so the retry actually fires instead of hanging on a dead
# mirror; and the loop falls through to named mirrors rather than trusting whichever host
# the ftpmirror redirector picks.
GSL_TARBALL="gsl-${GSL_VERSION}.tar.gz"
# Optional integrity check: set GSL_SHA256 to the published checksum to enable it.
GSL_SHA256="${GSL_SHA256:-}"

fetch_gsl() {
    curl --fail --location \
        --retry 5 --retry-delay 5 --retry-all-errors \
        --connect-timeout 30 --speed-limit 1024 --speed-time 30 \
        "$1" -o "${GSL_TARBALL}"
}

mkdir -p /tmp/gsl-src
cd /tmp/gsl-src
rm -f "${GSL_TARBALL}"
for url in \
    "https://ftpmirror.gnu.org/gsl/${GSL_TARBALL}" \
    "https://ftp.gnu.org/gnu/gsl/${GSL_TARBALL}" \
    "https://mirrors.kernel.org/gnu/gsl/${GSL_TARBALL}"
do
    echo "Fetching ${url}"
    if fetch_gsl "${url}"; then
        break
    fi
    echo "  mirror failed, trying next" >&2
    rm -f "${GSL_TARBALL}"
done

if [ ! -s "${GSL_TARBALL}" ]; then
    echo "ERROR: could not download ${GSL_TARBALL} from any mirror." >&2
    exit 1
fi

if [ -n "${GSL_SHA256}" ]; then
    echo "Verifying SHA-256..."
    if command -v sha256sum >/dev/null; then
        echo "${GSL_SHA256}  ${GSL_TARBALL}" | sha256sum --check --strict
    else
        echo "${GSL_SHA256}  ${GSL_TARBALL}" | shasum -a 256 --check --strict
    fi
else
    echo "GSL_SHA256 not set; skipping integrity check." >&2
fi

tar -xzf "gsl-${GSL_VERSION}.tar.gz"
cd "gsl-${GSL_VERSION}"

# Build & install
echo "Configuring GSL ${GSL_VERSION}..."
./configure --prefix="${INSTALL_DIR}" > /dev/null

echo "Building with ${CORES} cores..."
make -j"${CORES}" > /dev/null

echo "Installing..."
make install > /dev/null

echo "GSL ${GSL_VERSION} installed to ${INSTALL_DIR}"
