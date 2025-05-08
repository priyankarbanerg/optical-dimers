# Optical Dimers Simulation Codebase

This repository contains two main components for simulating optical dimers: the `cluster_version` and the `local_version`. The `cluster_version` is designed for running simulations on a computational cluster, while the `local_version` is used for running simulations locally on a single machine. Below are the details on how to use each version.

---

## **1. Cluster Version**

The `cluster_version` is designed to automate and run simulations sequentially on a computational cluster. It is well-integrated and handles tasks such as generating configurations, running simulations, applying convolutions, and saving results.

### **Key Features**
- Automates the simulation pipeline.
- Handles multiple configurations and instrument responses.
- Saves results in structured `.npy` files for further analysis.
- Designed for parallel execution on a cluster.

### **How to Use**
1. **Set Up Configurations**:
   - Configurations are defined in the `run.py` or `run_old.py` files.
   - Example configurations:
     ```python
     configurations = [
         {"config_name": "H", "faulty_detector": "no", "tumbling": "slow", "jiggling": "no", "mode": "polaron"},
         {"config_name": "J", "faulty_detector": "yes", "tumbling": "fast", "jiggling": "yes", "mode": "weak"}
     ]
     ```

2. **Run Simulations**:
   - Use the `G2InstrumentResponseSimulator` or `G2CorrelationSimulator` classes to run simulations.
   - Example:
     ```python
     instrument_simulator = G2InstrumentResponseSimulator(configurations=configurations)
     instrument_simulator.run_simulation_for_responses()
     ```

3. **Output**:
   - Results are saved in the `files/instrument_response` or `files/g2_con` directories.
   - File naming conventions include configuration details (e.g., `H_strongC_2nm_symm_slow_tumb_jigg_convolved.npy`).

4. **Run Commands**:
   - Submit the script to the cluster using your cluster's job submission system (e.g., `sbatch` for SLURM).

---

## **2. Local Version**

The `local_version` is designed for running simulations on a local machine. It is less automated and requires manual intervention for certain tasks such as plotting and data analysis.

### **Key Features**
- Allows for quick testing and visualization of results.
- Includes plotting utilities for analyzing simulation outputs.
- Requires manual setup and execution of tasks.

### **How to Use**
1. **Run Simulations**:
   - Use the `propagate.py` file to run simulations.
   - Example:
     ```python
     # Run a specific simulation
     python local_version/propagate.py
     ```

2. **Plot Results**:
   - Use the `plots.py` file to visualize results.
   - Example:
     ```python
     python local_version/plots.py
     ```

3. **Output**:
   - Results are saved in the `files/g2` or `files/spectra` directories.
   - File naming conventions are similar to the cluster version.

4. **Manual Adjustments**:
   - Modify parameters directly in the scripts (e.g., file paths, configurations, plotting ranges).

---

## **File Structure**
- **`cluster_version/`**:
  - `run.py`: Main script for running simulations on the cluster.
  - `run_old.py`: Legacy script with additional configurations and methods.
  - `system_ops.py`: Contains utility functions for system operations.
  - `files/`: Directory for saving simulation outputs.

- **`local_version/`**:
  - `propagate.py`: Script for running local simulations.
  - `plots.py`: Script for visualizing results.
  - `files/`: Directory for saving local outputs.

---

## **General Notes**
- Ensure all dependencies are installed before running the scripts (e.g., `numpy`, `matplotlib`).
- For cluster execution, ensure the cluster environment supports Python and the required libraries.
- Modify configurations and parameters as needed for your specific use case.

---

## **Future Improvements**
- Integrate the `local_version` with the `cluster_version` for better automation.
- Add detailed logging and error handling for both versions.
- Provide a unified interface for running simulations locally and on the cluster.

---

## **Contact**
For any issues or questions, please contact the repository owner.
