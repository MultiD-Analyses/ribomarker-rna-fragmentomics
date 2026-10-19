# Licence texts added to the image

The conda packages below do not include their licence text, so the texts are taken from the upstream source at the packaged version and copied into the image (`/opt/conda/envs/ribomarker/share/conda-licenses/_added/`) by the Dockerfile. Downloaded 2026-10-15.

| File | Package(s) | Source |
|---|---|---|
| `bowtie-1.0.0_LICENSE` | bowtie 1.0.0 | https://github.com/BenLangmead/bowtie/blob/v1.0.0/LICENSE (the Artistic License; the conda package records GPLv3) |
| `bedtools-2.31.1_LICENSE` | bedtools 2.31.1 | https://github.com/arq5x/bedtools2/blob/v2.31.1/LICENSE |
| `dnaio-1.2.3_LICENSE` | dnaio 1.2.3 | https://github.com/marcelm/dnaio/blob/v1.2.3/LICENSE |
| `pysam-0.24.0_COPYING` | pysam 0.24.0 | https://github.com/pysam-developers/pysam/blob/v0.24.0/COPYING |
| `gcc-15.2.0_COPYING3`, `gcc-15.2.0_COPYING.RUNTIME` | libgcc, libgcc-ng, libgfortran, libgfortran5, libgomp, libstdcxx, libstdcxx-ng 15.2.0 | https://github.com/gcc-mirror/gcc/tree/releases/gcc-15.2.0 (GPL v3 with the GCC Runtime Library Exception 3.1) |

Not added: libfreetype and libfreetype6 2.14.3 (their texts are in the `freetype-2.14.3-*` folder of the same source), libsqlite and sqlite 3.53.1 (public domain).
