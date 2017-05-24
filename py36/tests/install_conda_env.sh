
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
export PATH="~/miniconda2/bin:$PATH"
source ~/.bashrc
conda --version
echo "channels:" > ~/.condarc
echo "  - bioconda" >> ~/.condarc
echo "  - r" >> ~/.condarc
echo "  - defaults" >> ~/.condarc
. activate RepEnrich_py36_0

