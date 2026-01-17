# ECDH — Elliptic Curve Diffie–Hellman (C++)

This repository contains an educational and research-oriented implementation
of the Elliptic Curve Diffie–Hellman (ECDH) key exchange protocol written in C++.

The project is developed from scratch without using third-party cryptographic
libraries (OpenSSL, libsodium, etc.) and focuses on understanding the internal
mechanisms of elliptic curve cryptography.

The work is accompanied by a scientific paper and analytical experiments and
was presented at the international conference:

**Optical Technologies, Materials & Systems 2025**  
RTU MIREA  
https://www.mirea.ru/conference-optical-technologies-materials-systems/2025/

Repository:  
https://github.com/GeeZiX84/ECDH

---

## Documentation

- [Scientific paper (interactive)](docs/paper.md)
- [Theory vs Implementation](docs/theory_vs_code.md)
- [Differences: Theory vs Code](docs/differences_theory_code.md)
- [Analytical report (ECDH vs BD, PDF)](docs/analysis_ecdh_bd.pdf)
- [Python analysis](docs/analysis_ecdh_bd.pdf)
- [README на русском](README_RU.md)

---


## Project Goals

- implement finite field arithmetic over prime fields
- implement elliptic curve group operations
- demonstrate the ECDH key exchange protocol
- theoretically and experimentally compare ECDH and Burmester–Desmedt protocols
- connect cryptographic theory with C++ and Python implementations

---

## Repository Structure

ECDH/
├── include/ # C++ headers
├── src/ # Core C++ implementation
├── python/ # Python analytical code
├── docs/ # Scientific documentation
├── examples/ # Optional demo programs
└── CMakeLists.txt

---

## Documentation

- [Interactive scientific paper](docs/paper.md)
- [Analytical report: ECDH vs BD (PDF)](docs/analysis_ecdh_bd.pdf)
- [Conference presentation (PDF)](docs/presentation.pdf)
- [Theory vs Implementation](docs/theory_vs_code.md)

---

## Build Requirements

### C++ Implementation

- C++17 (explicitly tested with GCC via MinGW-w64)
- CMake >= 3.16

Supported environments:
- Linux (GCC / Clang)
- Windows x64 (MinGW-w64 with Windows SDK headers)

The project relies on Windows SDK headers provided by MinGW-w64
(e.g. `ncrypt.h`, `bcrypt.h`, `cert*`, `cfgmgr32.h`) when built on Windows and headers from [boost/multiprecision](https://github.com/boostorg/multiprecision), [nlohmann/json](https://github.com/J08nY/ecgen?tab=readme-ov-file). JSON by nlohman was used for reading json that was generated with [J08nY's "ecgen"](https://github.com/J08nY/ecgen?tab=readme-ov-file).

No third-party cryptographic libraries are used.

---

## Build Instructions

### Linux

```bash
git clone https://github.com/GeeZiX84/ECDH.git
cd ECDH

cmake -S . -B build
cmake --build build
```

Windows (MinGW-w64)

Expected MinGW-w64 installation path: "C:/msys64/mingw64"
```bash
git clone https://github.com/GeeZiX84/ECDH.git
cd ECDH

cmake -G "MinGW Makefiles" -S . -B build
cmake --build build
```

---

## Notes and Limitations

1. arithmetic is implemented explicitly

2. no cryptographically secure RNG is provided

3. no side-channel attack countermeasures

4. Python code is used only for analytical evaluation

This project must not be used in production cryptographic systems.

[In realisation of this project was used tool "ecgen" for generating ecliptic curves by J08nY](https://github.com/J08nY/ecgen?tab=readme-ov-file')

---

## License

Distributed for educational and research purposes.


---