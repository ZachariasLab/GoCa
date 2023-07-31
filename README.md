# Generate structure-based coarse-grained protein models for molecular dynamics simulations

### Installation
```bash
git clone git@gitlab.lrz.de:luiswalter/coarse-grained-protein-simulation.git
cd coarse-grained-protein-simulation
git submodule init # Fetch dependencies
git submodule update
make
```

### Usage
1. Copy the default configuration file `default.ini` to your target directory and modify it according to your needs.
2. Run the following command to generate a topology and coordinate file. If no configuration file is provided the script looks for a file called `GoCa.ini`.
   ```
   ./path_to_repository/GoCa GoCa.ini
   ```
3. Use the topology and coordinate files for your simulation

### Tests
To run the provided test install `pytest` and run from the main directory
```
python3 -m pytest tests
```