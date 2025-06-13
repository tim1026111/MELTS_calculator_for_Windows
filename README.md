# MELTS Calculator - Windows Edition

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A Python-based GUI application for thermodynamic equilibrium calculations in magmatic systems, based on the MELTS algorithm. This Windows-compatible version provides an intuitive interface for petrology and geochemistry calculations.

## 🌟 Features

- **User-friendly GUI** with tabbed interface
- **Thermodynamic calculations** for magmatic systems
- **Equilibrium and fractional crystallization** modeling
- **Real-time visualization** of results
- **Multiple phase saturation** calculations
- **Temperature-pressure sequence** computations
- **Data import/export** functionality
- **Offline operation** - no internet required

## 📋 Requirements

- Windows 10 or higher
- Python 3.8 or higher
- 4GB RAM (minimum)
- 100MB disk space

## 🚀 Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/tim1026111/melts-calculator.git
cd melts-calculator
```

### 2. Install dependencies
```bash
pip install -r requirements.txt
```

### 3. Run the application
```bash
python melts_calculator.py
```

## 📦 Installation

### Option 1: Using pip
```bash
# Install required packages
pip install numpy scipy pandas matplotlib

# Download the repository
git clone https://github.com/tim1026111/melts-calculator.git
cd melts-calculator

# Run the program
python melts_calculator.py
```

### Option 2: Using Anaconda
```bash
# Create a new environment
conda create -n melts python=3.8
conda activate melts

# Install packages
conda install numpy scipy pandas matplotlib

# Run the program
python melts_calculator.py
```

## 📖 Usage Guide

### 1. Input Composition
- Navigate to the **Composition** tab
- Enter oxide weight percentages
- Click "Normalize to 100%" to normalize values
- Molar values are calculated automatically

### 2. Set Conditions
- Go to the **Conditions** tab
- Set temperature (°C), pressure (bar), and oxygen fugacity
- Choose calculation mode (equilibrium or fractional)
- Optional: Enable T-P sequence for path calculations

### 3. Select Phases
- In the **Phases** tab, select minerals to include
- Default phases: liquid, feldspar, quartz
- Use "Select All" or "Deselect All" for convenience

### 4. Run Calculations
- Click the **Calculate** button
- View results in the **Results** tab
- Check visualizations in the **Plots** tab

### 5. Save/Load Data
- Save results as JSON files
- Load previous calculations
- Export graphs as images

## 📊 Example Calculation

### Granite Composition
```
SiO2:  74.39 wt%
Al2O3: 13.55 wt%
FeO:   0.98 wt%
MgO:   0.50 wt%
CaO:   1.43 wt%
Na2O:  3.36 wt%
K2O:   5.09 wt%
H2O:   5.04 wt%
```

### Conditions
- Temperature: 800°C
- Pressure: 200 bar
- log fO2: 0 (QFM buffer)

## 🔧 Troubleshooting

### Common Issues

1. **"Python is not recognized"**
   - Ensure Python is added to PATH
   - Restart command prompt after Python installation

2. **Import errors**
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt
   ```

3. **GUI not displaying properly**
   - Check display scaling settings
   - Try running in compatibility mode

## 🏗️ Project Structure

```
melts-calculator/
│
├── melts_calculator.py      # Main GUI application
├── melts_thermodynamics.py  # Thermodynamic calculations
├── requirements.txt         # Python dependencies
├── README.md               # This file
├── LICENSE                 # MIT License
├── examples/               # Example data files
│   └── granite_example.json
└── docs/                   # Additional documentation
    └── theory.md          # Theoretical background
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📚 Background

This calculator is based on the MELTS algorithm developed by Ghiorso and Sack (1995). The implementation provides a simplified version suitable for educational and basic research purposes.

### Key References
- Ghiorso, M.S. and Sack, R.O. (1995) Chemical Mass Transfer in Magmatic Processes. IV. A Revised and Internally Consistent Thermodynamic Model for the Interpolation and Extrapolation of Liquid-Solid Equilibria in Magmatic Systems at Elevated Temperatures and Pressures. Contributions to Mineralogy and Petrology, 119, 197-212.

## ⚠️ Disclaimer

This is a simplified implementation of the MELTS algorithm. For research-grade calculations, please use the original MELTS software or consult with experts in the field.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- Byeon jeong bin - Initial work

## 🙏 Acknowledgments

- Original MELTS algorithm by Mark Ghiorso and Richard Sack
- Python scientific computing community
- Contributors and testers
