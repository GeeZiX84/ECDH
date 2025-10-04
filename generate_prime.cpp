#include <iostream>
#include <fstream>
#include <random>
#include <cstdint>

// проверка на простоту простым перебором (для малых чисел хватит)
bool is_prime(int64_t n) {
    if (n < 2) return false;
    for (int64_t i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

int main() {
    std::ofstream fout("primes.txt");
    if (!fout) {
        std::cerr << "Ошибка: не удалось открыть файл\n";
        return 1;
    }

    int64_t limit = 1000000; // до какого числа генерировать
    for (int64_t i = 1; i <= limit; i++) {
        if (is_prime(i)) {
            fout << i << "\n"; // каждое простое в новой строке
        }
    }

    std::cout << "Простые числа сохранены в primes.txt\n";
}
