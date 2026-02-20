# syntax=docker/dockerfile:1
FROM mcr.microsoft.com/devcontainers/anaconda:1-3

# Copy environment file and install conda environment
COPY leafcutter2_env.yml /tmp/leafcutter2_env.yml
RUN conda env create -f /tmp/leafcutter2_env.yml && rm /tmp/leafcutter2_env.yml

# Set conda environment as default
SHELL ["/bin/bash", "-c"]
ENV CONDA_DEFAULT_ENV=leafcutter2
ENV PATH=/opt/conda/envs/leafcutter2/bin:$PATH

# Activate environment for all bash sessions
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate leafcutter2" >> /etc/bash.bashrc

# Set workdir
WORKDIR /workspace

COPY . /workspace

# Set the default command to run leafcutter2.py, passing any arguments
ENTRYPOINT ["python", "/workspace/scripts/leafcutter2.py"]
# To pass parameters, use: docker run <image> -j junction_files.txt ...
