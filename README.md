# ECG QRS Detection Methods

This repository contains MATLAB implementations of three different methods for QRS complex detection in ECG signals. The project uses the MIT-BIH Arrhythmia Database for testing and evaluation.

## Methods Implemented

1. Zero-Crossing Method
2. Empirical Mode Decomposition (EMD) Method
3. Morphological Operations Method

## Requirements

- MATLAB (version R2019b or later recommended)
- [WFDB Toolbox for MATLAB](https://physionet.org/content/wfdb-matlab/0.10.0/)
- [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/)

## File Structure

### Main Scripts
- `main_crossing.m`: Main script for QRS detection using the zero-crossing method
- `main_emd.m`: Main script for QRS detection using the EMD method
- `main_morphological.m`: Main script for QRS detection using morphological operations

### Detection Functions
- `detect_qrs_crossing.m`: Implements the zero-crossing method for QRS detection
- `detect_qrs_emd.m`: Implements the EMD method for QRS detection
- `detect_qrs_morphological.m`: Implements the morphological operations method for QRS detection

### Utility Functions
- `evaluate_detection.m`: Function to evaluate the performance of QRS detection algorithms

## Setup

1. Clone this repository to your local machine.
2. Download the [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/) and place it in a folder named `mit-bih-arrhythmia-database-1.0.0` in the project directory.
3. Install the [WFDB Toolbox for MATLAB](https://physionet.org/content/wfdb-matlab/0.10.0/) and ensure it's in your MATLAB path.

## Usage

1. Open MATLAB and navigate to the project directory.
2. Run one of the main scripts:
   - `main_crossing.m`
   - `main_emd.m`
   - `main_morphological.m`
3. The script will process the ECG signal, detect QRS complexes, and display the results along with performance metrics (if annotations are available).

## Customization

- To analyze a different record from the MIT-BIH Arrhythmia Database, change the `record_name` variable in the main script you're using.

## Performance Metrics

The scripts calculate and display two performance metrics when annotations are available:

1. Sensitivity (Se): The proportion of actual QRS complexes that were correctly detected.
2. Positive Predictivity (+P): The proportion of detected QRS complexes that were actual QRS complexes.

These metrics are calculated using the `evaluate_detection.m` function.

## Contributing

Contributions to improve the algorithms or add new QRS detection methods are welcome. Please feel free to submit pull requests or open issues for any bugs or suggestions.

## License

This project is open-source and available under the MIT License.

## Acknowledgments

- This project uses the [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/) from PhysioNet.
- The WFDB Toolbox for MATLAB is used for reading ECG signals and annotations.
