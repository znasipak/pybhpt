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
# curl's own --retry only covers transient HTTP responses (408, 429, 5xx) and timeouts, so
# a mid-transfer connection reset -- exit 35, which is how a bad GNU mirror usually fails --
# is not retried. --retry-all-errors would fix that but needs curl >= 7.71, and the
# manylinux/musllinux build containers ship far older curl (CentOS 7 is on 7.29), so the
# retrying is done here in shell instead: every attempt is retried across every mirror.
# Only long-standing curl options are used. --speed-limit/--speed-time abort a transfer
# that has stalled so the next attempt starts instead of hanging on a dead mirror.
GSL_TARBALL="gsl-${GSL_VERSION}.tar.gz"
# Optional integrity check: set GSL_SHA256 to the published checksum to enable it.
GSL_SHA256="${GSL_SHA256:-}"

GSL_MIRRORS="
https://ftpmirror.gnu.org/gsl/${GSL_TARBALL}
https://ftp.gnu.org/gnu/gsl/${GSL_TARBALL}
https://mirrors.kernel.org/gnu/gsl/${GSL_TARBALL}
"
GSL_FETCH_ATTEMPTS="${GSL_FETCH_ATTEMPTS:-3}"

fetch_gsl() {
    curl --fail --location \
        --connect-timeout 30 --speed-limit 1024 --speed-time 30 \
        "$1" -o "${GSL_TARBALL}"
}

mkdir -p /tmp/gsl-src
cd /tmp/gsl-src
rm -f "${GSL_TARBALL}"
attempt=1
while [ "${attempt}" -le "${GSL_FETCH_ATTEMPTS}" ]; do
    for url in ${GSL_MIRRORS}; do
        echo "Fetching ${url} (attempt ${attempt}/${GSL_FETCH_ATTEMPTS})"
        if fetch_gsl "${url}"; then
            break 2
        fi
        echo "  mirror failed" >&2
        rm -f "${GSL_TARBALL}"
    done
    attempt=$((attempt + 1))
    if [ "${attempt}" -le "${GSL_FETCH_ATTEMPTS}" ]; then
        echo "  all mirrors failed; retrying in 5s" >&2
        sleep 5
    fi
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
