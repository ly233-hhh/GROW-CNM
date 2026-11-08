# GROW-CNM
Source code for Global weight-Ranked cluster-grOWing (GROW) algorithm and CNM pipeline.

# Global weight-Ranked cluster-grOWing (GROW)

This repository contains the official MATLAB implementation of the **Global weight-Ranked cluster-grOWing (GROW)** algorithm and the associated **Coordinate Network Mapping (CNM)** pipeline. 

GROW is designed to reconstruct the spatial extent of statistical clusters from reported peak coordinates and a reference statistical weight map, without requiring original raw imaging data. It replaces the conventional fixed-radius spherical approximation in CNM, providing better spatial fidelity and robustness to study-level confounders.

## 📋 Prerequisites

To run these scripts, you must have:
*   **MATLAB** (Tested on R2020a and newer versions)
*   **SPM12** (Statistical Parametric Mapping, for reading/writing NIfTI files: `spm_vol`, `spm_read_vols`, `spm_write_vol`). Please ensure SPM12 is added to your MATLAB path.

## 📁 File Structure

*   `generate_cluster.m`: The core algorithm script for reconstructing a single cluster based on a seed MNI coordinate, a target cluster size, and a weight map (Best-First Region Growing).
*   `bqw_networkmapping_weight2.m`: The main network mapping pipeline that processes multiple seeds, generates virtual clusters, and computes functional connectivity (FC) maps using a normative connectome.
*   `coord_mni2vox.m`: Utility function for coordinate transformation between MNI space and voxel space.
*   `SphereMask.m`: Utility function for generating traditional fixed-radius spherical masks (used for comparison with Fixed-CNM).
*   `virtual_cluster_gen_test.m`: An example script demonstrating how to use `generate_cluster.m` with a mock/example weight map.

## 🚀 Quick Start

1.  Clone this repository or download the `.m` files.
2.  Open MATLAB and navigate to the downloaded folder.
3.  Add the folder to your MATLAB path: `addpath(genpath(pwd));`
4.  Add SPM12 to your path.
5.  Open `virtual_cluster_gen_test.m`, modify the `weight_map_path` and `mni_coord` to point to your own data, and run the script.

## 📖 Citation

If you use GROW or this code in your research, please cite our paper:

> Liu Y, Fu X, Wang X, et al. Global weight-Ranked cluster-grOWing (GROW): A Cluster-Growing Framework for Robust Functional Coordinate Network Mapping. *NeuroImage* (Under Review).

*(Citation details will be updated upon publication).*

## 📧 Contact

For questions regarding the code or the algorithm, please contact:
*   **Wen Qin** - [wayne.wenqin@gmail.com]
*   **Yu Liu** - [17347815310@163.com]

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
