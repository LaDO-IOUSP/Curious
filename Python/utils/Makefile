run: conda-packages, OceanLab, dedalus

# make conda-packages
conda-packages:
	yes | conda install -c conda-forge pandas
	yes | conda install -c conda-forge seawater
	yes | conda install -c conda-forge gsw
	yes | conda install -c conda-forge netcdf4
	yes | conda install -c conda-forge pyproj
	yes | conda install -c conda-forge proj4
	yes | conda install -c conda-forge basemap-data-hires
	yes | conda install -c conda-forge basemap
	yes | conda install -c conda-forge cmocean
	yes | conda install -c conda-forge seaborn

# make OceaLab AT=<installation directory>
OceanLab:
	# choose installation path
	@echo OceanLab installation at $(AT)
	if [ ! -d $$AT ]; then mkdir $$AT; fi
	git clone https://github.com/iuryt/OceanLab.git $$AT
	cd $$AT && python setup.py install

# make dedalus AT=<installation directory>
dedalus:
	yes | conda install -n base -c conda-forge mpi
	yes | conda install -n base -c conda-forge mpi4py
	yes | conda install -n base -c conda-forge openmpi
	yes | conda install -n base -c conda-forge fftw
	yes | conda install -n base -c conda-forge hdf5
	# choose installation path
	@echo dedalus installation at $(AT)
	if [ ! -d $$AT ]; then mkdir $$AT; fi
	# download dedalus installer
	wget https://raw.githubusercontent.com/DedalusProject/conda_dedalus/master/install_conda.sh -O $$AT/install_conda.sh
	# install dedalus	
	bash $$AT/install_conda.sh
	# activate dedalus and install packages
	reset
	yes | conda install -n dedalus -c anaconda ipython
	yes | conda install -n dedalus -c conda-forge seawater
	yes | conda install -n dedalus -c conda-forge pandas
	yes | conda install -n dedalus -c conda-forge seaborn

# make parcels ENV=<env_name>
parcels:
	yes | conda create -n $$ENV -c conda-forge python=3.6 parcels jupyter cartopy ffmpeg
