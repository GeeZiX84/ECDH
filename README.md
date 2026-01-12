# ECDH — Elliptic Curve Diffie–Hellman (C++)

This repository contains an educational and research-oriented implementation
of the Elliptic Curve Diffie–Hellman (ECDH) key exchange protocol written in C++.

The project is developed from scratch without using third-party cryptographic
libraries (OpenSSL, libsodium, etc.) and is intended to demonstrate the internal
mechanisms of elliptic curve cryptography.

The implementation is accompanied by a scientific paper presented at the
conference:

Optical Technologies, Materials & Systems 2025  
RTU MIREA  
https://www.mirea.ru/conference-optical-technologies-materials-systems/2025/

GitHub repository:  
https://github.com/GeeZiX84/ECDH

---

## Project Goals

- implement finite field arithmetic over prime fields
- implement elliptic curve group operations
- demonstrate the ECDH key exchange protocol
- connect theoretical cryptographic analysis with practical implementation
- provide a clean and portable CMake-based build system

---

## Project Structure

ECDH/
├── include/ # Public headers
├── src/ # Core implementation
├── examples/ # Optional demo applications
├── docs/ # Paper, presentation and documentation
├── CMakeLists.txt
└── README.md


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
Windows (MinGW-w64)
Expected MinGW-w64 installation path:


C:/msys64/mingw64
Build commands:


git clone https://github.com/GeeZiX84/ECDH.git
cd ECDH

cmake -G "MinGW Makefiles" -S . -B build
cmake --build build
```

System headers from the MinGW-w64 Windows SDK
(ncrypt.h, cert*, cfgmgr32.h, etc.) are automatically detected by CMake.

---

## Implementation Notes
1. All arithmetic is implemented explicitly
2. No platform-independent cryptographic RNG is provided 
3. No countermeasures against side-channel attacks are implemented
4. Windows-specific headers may be used when building under MinGW-w64

This project must not be used in production cryptographic systems.

---
## Scientific Context
The implementation is based on a theoretical analysis of key exchange algorithms over elliptic curves. In addition to ECDH, the Burmester–Desmedt group key exchange protocol is
analyzed theoretically for comparison purposes.

---

## License
The project is distributed for educational and research purposes.

---

## ECDH — протокол Диффи–Хеллмана на эллиптических кривых (C++)

Данный репозиторий содержит учебную и исследовательскую реализацию протокола
обмена ключами Elliptic Curve Diffie–Hellman (ECDH), написанную на языке C++.

Проект реализован «с нуля» без использования сторонних криптографических
библиотек (OpenSSL, libsodium и др.) и предназначен для изучения внутренних
механизмов эллиптической криптографии.

Работа сопровождается научной статьёй и была представлена на конференции:

**Optical Technologies, Materials & Systems 2025**  
РТУ МИРЭА  
https://www.mirea.ru/conference-optical-technologies-materials-systems/2025/

Репозиторий проекта:  
https://github.com/GeeZiX84/ECDH

---

## Цели проекта

- реализация арифметики конечных полей над простым модулем  
- реализация групповых операций на эллиптических кривых  
- демонстрация протокола обмена ключами ECDH  
- связь теоретического криптографического анализа с практической реализацией  
- использование корректной CMake-сборки для разных платформ  

---

## Структура проекта

ECDH/
├── include/ 
├── src/ 
├── examples/  
├── docs/ 
├── CMakeLists.txt
└── README.md

---

## Требования для сборки

- стандарт C++17 или новее  
- CMake версии 3.16 или выше  
- компилятор:
  - Linux: GCC или Clang  
  - Windows: MinGW-w64 (x86_64)

---

## Сборка проекта

### Linux

```bash
git clone https://github.com/GeeZiX84/ECDH.git
cd ECDH

cmake -S . -B build
cmake --build build
Windows (MinGW-w64)
Ожидаемый путь установки MinGW-w64:

C:/msys64/mingw64
Сборка выполняется следующими командами:

git clone https://github.com/GeeZiX84/ECDH.git
cd ECDH

cmake -G "MinGW Makefiles" -S . -B build
cmake --build build
```

При сборке под Windows автоматически используются системные заголовочные файлы
MinGW-w64 (включая ncrypt.h, cert*, cfgmgr32.h и др.).

---

## Особенности и ограничения реализации
1. Арифметика конечных полей реализована вручную

2. Отсутствует криптографически стойкий генератор случайных чисел

3. Не реализованы защиты от атак по побочным каналам

4. Возможна платформозависимость при использовании Windows SDK

Проект не предназначен для использования в промышленной криптографии и рассматривается исключительно как учебный и исследовательский.

---

## Научный контекст
В рамках работы проведён теоретический анализ алгоритмов обмена ключами на эллиптических кривых. Помимо ECDH, в статье рассматривается протокол группового обмена ключами Burmester–Desmedt (теоретически, без программной реализации).

---

## Лицензия
Проект распространяется в учебных и научных целях.

---

