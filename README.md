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

- C++17 or newer
- CMake >= 3.16
- Compiler:
  - Linux: GCC or Clang
  - Windows: MinGW-w64 (x86_64)

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

---

## License

Distributed for educational and research purposes.


---