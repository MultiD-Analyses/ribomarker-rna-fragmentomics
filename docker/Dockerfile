# Image of the pipeline (publication plan, step 2): one stage, plain Python scripts, the exact conda
# environment (env.lock.txt, with builds) and the human and mouse references baked in.
#
# Build (from the repository root; references staged first by docker/stage_references.py):
#   docker build --provenance=false \
#     --build-context references=staging/references --build-context root=. \
#     --build-arg GIT_COMMIT=$(git rev-parse HEAD) --build-arg VERSION=$(cat docker/VERSION) \
#     -t IMAGE:TAG docker
#
# Base image: Miniforge (conda-forge only), pinned by digest (condaforge/miniforge3:latest, pulled 2026-10-13).
FROM condaforge/miniforge3@sha256:f752860f77bd417aa4db35be3d2e906fdb740fb80f40755ed651fff9e8873aad

WORKDIR /pipeline

# References first: largest and rarely changed, so a code change does not rebuild these layers.
# Only the files listed in MANIFEST_human.json and MANIFEST_mouse.json (plus MANIFESTs and PROVENANCEs), checked against their
# checksums by the staging script. Each large genome index file gets its own layer (upload limits).
COPY --from=references rest/ /pipeline/assets/
COPY --from=references big/genome/human-genome.1.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/human-genome.2.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/human-genome.3.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/human-genome.4.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/human-genome.rev.1.bt2 /pipeline/assets/genome/
COPY --from=references big/genome/human-genome.rev.2.bt2 /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.1.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.2.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.3.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.4.bt2     /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.rev.1.bt2 /pipeline/assets/genome/
COPY --from=references big/genome/mouse-genome.rev.2.bt2 /pipeline/assets/genome/

# The conda environment: exactly the validated package builds. conda clean in the same step keeps
# the package cache out of the layer. The licence texts live in that cache (pkgs/<package>/info/licenses),
# so they are copied first to /opt/conda/envs/ribomarker/share/conda-licenses/<package>/ (GPL needs them).
COPY env.lock.txt /pipeline/env.lock.txt
RUN conda create -y -n ribomarker --file /pipeline/env.lock.txt \
 && mkdir -p /opt/conda/envs/ribomarker/share/conda-licenses \
 && for d in /opt/conda/pkgs/*/info/licenses; do \
      p=${d%/info/licenses}; cp -r "$d" "/opt/conda/envs/ribomarker/share/conda-licenses/${p##*/}"; \
    done \
 && conda clean -afy
# Licence texts of the few packages that do not include their own (sources: licenses/README.md).
COPY licenses/ /opt/conda/envs/ribomarker/share/conda-licenses/_added/

# The references' <species>.dbpaths.csv holds paths relative to its folder, so the same file works
# here and locally. RIBOMARKER_CONTAINER=1 makes the pipeline refuse older absolute-path dbpaths files.
ENV RIBOMARKER_CONTAINER=1
ENV RIBOMARKER_PIPELINE_DIR=/pipeline

# Demo data for `run --demo`: small synthetic reads, truth and references (made by docker/make_demo.sh).
COPY demo/ /pipeline/demo/

# Pipeline code (changes often): plain Python scripts, the Snakefile, QC settings, the RNA class tables
# (refTables.py), the entry script.
COPY scripts/ /pipeline/scripts/
COPY qc/ /pipeline/qc/
COPY tables/ /pipeline/tables/
COPY Snakefile VERSION run-pipeline.sh /pipeline/
RUN chmod +x /pipeline/run-pipeline.sh \
 && /opt/conda/envs/ribomarker/bin/python /pipeline/scripts/codeFingerprint.py /pipeline > /pipeline/CODE_FINGERPRINT

# Licence, notices and third-party notices travel with the image (MIT notice; data attributions).
COPY --from=root LICENSE NOTICE THIRD_PARTY_NOTICES.md /pipeline/

# Version and commit of this build (shown in the QC report and as image labels).
ARG GIT_COMMIT=unknown
ARG VERSION=unknown
RUN echo "$GIT_COMMIT" > /pipeline/GIT_COMMIT
LABEL org.opencontainers.image.source="https://github.com/MultiD-Analyses/ribomarker-rna-fragmentomics" \
      org.opencontainers.image.version="$VERSION" \
      org.opencontainers.image.revision="$GIT_COMMIT"

ENTRYPOINT ["/pipeline/run-pipeline.sh"]
