#!/usr/bin/env bash
set -euo pipefail

EX47_SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
EX47_BUILD_DIR="${EX47_BUILD_DIR:-${EX47_SOURCE_DIR}/build}"
EX47_CMAKE="${CMAKE:-cmake}"
EX47_CONFIGURE=(-S "${EX47_SOURCE_DIR}" -B "${EX47_BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release)
if [[ -n "${EX47_BOOST_INCLUDE_DIR:-}" ]]; then
  EX47_CONFIGURE+=("-DEX47_BOOST_INCLUDE_DIR=${EX47_BOOST_INCLUDE_DIR}")
fi
if [[ -n "${EX47_EIGEN_INCLUDE_DIR:-}" ]]; then
  EX47_CONFIGURE+=("-DEX47_EIGEN_INCLUDE_DIR=${EX47_EIGEN_INCLUDE_DIR}")
fi
"${EX47_CMAKE}" "${EX47_CONFIGURE[@]}"
"${EX47_CMAKE}" --build "${EX47_BUILD_DIR}" --parallel "${EX47_JOBS:-2}"
exec "${EX47_BUILD_DIR}/MultiPhase_ex47" "$@"
