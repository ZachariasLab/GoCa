# GoCa
**Generate structure-based coarse-grained protein models for molecular dynamics simulations!**  
The GoCa program can also be used via the [luiswalter.dev/GoCa](https://luiswalter.dev/GoCa/) webserver.  

### Installation
We use Make to build the binary.
```bash
git clone git@github.com:daandtu/GoCa.git
cd GoCa
git submodule init # Fetch dependencies
git submodule update
make
```

### Usage
1. Copy the default configuration file `default.ini` to your target directory and modify it according to your needs.
2. Execute the program binary to generate a topology and coordinate file. Pass the previously created configuration file as a command line argument. If no configuration file is specified, the script searches for a file named `GoCa.ini`.
   ```
   ./path_to_this_repository/GoCa configuration.ini
   ```
3. Use the topology and coordinate files for your simulation (for example, with [this](https://luiswalter.dev/GoCa/gromacs-example) GROMACS configuration).

### Tests
To run the provided test, install `pytest` and run from the main directory
```
python3 -m pytest tests
```
