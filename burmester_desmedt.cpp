// bd_example.cpp
// Демонстрация: поиск генератора g и реализация протокола Burmester-Desmedt (BD).
// Компиляция: g++ -std=c++17 -O2 bd_example.cpp -o bd_example

#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstdint>

using boost::multiprecision::cpp_int;
using std::vector;
using std::uint64_t;
using std::cout;
using std::endl;

// Возведение a^e mod m (быстрое возведение для cpp_int)
cpp_int mod_pow(cpp_int a, cpp_int e, const cpp_int &m) {
    cpp_int res = 1;
    a %= m;
    while (e > 0) {
        if ((e & 1) != 0) res = (res * a) % m;
        a = (a * a) % m;
        e >>= 1;
    }
    return res;
}

// Евклид - обратный элемент a^{-1} mod m (расширенный алгоритм Евклида).
// Возвращает (inv, gcd). Если gcd != 1 — обратного нет.
std::pair<cpp_int, cpp_int> egcd(cpp_int a, cpp_int b) {
    if (b == 0) return {1, a};
    cpp_int x0 = 1, x1 = 0;
    cpp_int y0 = 0, y1 = 1;
    while (b != 0) {
        cpp_int q = a / b;
        cpp_int r = a % b;
        cpp_int x2 = x0 - q * x1;
        cpp_int y2 = y0 - q * y1;
        a = b; b = r;
        x0 = x1; x1 = x2;
        y0 = y1; y1 = y2;
    }
    // a == gcd, coefficients x0,y0: x0*origA + y0*origB = gcd
    return {x0, a};
}

bool inv_mod(cpp_int a, const cpp_int &m, cpp_int &out_inv) {
    auto p = egcd(a, m);
    cpp_int inv = p.first;
    cpp_int g = p.second;
    if (g != 1 && g != -1) return false;
    inv %= m;
    if (inv < 0) inv += m;
    out_inv = inv;
    return true;
}

// Простое пробное деление (подходит для небольших p-1).
// Возвращает список простых множителей (без повторов).
vector<cpp_int> factorize_small(const cpp_int &n) {
    vector<cpp_int> res;
    cpp_int x = n;
    // делим на 2
    if (x % 2 == 0) {
        res.push_back(2);
        while (x % 2 == 0) x /= 2;
    }
    // пробное деление по нечетным до некоторого предела
    uint64_t limit = 2000000; // можно увеличить или уменьшить
    for (uint64_t d = 3; d <= limit && x > 1; d += 2) {
        if (x % d == 0) {
            res.push_back(d);
            while (x % d == 0) x /= d;
        }
    }
    if (x > 1) {
        // остался составной или большой простой множитель; добавим как "фактор"
        res.push_back(x);
    }
    return res;
}

// Поиск генератора g для группы порядка p (попытается факторизовать p-1).
// Возвращает 0, если не найден (маловероятно для разумных p).
cpp_int find_generator(const cpp_int &p) {
    if (p <= 3) return 0;
    cpp_int phi = p - 1;
    auto primes = factorize_small(phi);
    // попробуем кандидатов g = 2,3,4,...
    for (cpp_int g = 2; g < p; ++g) {
        bool ok = true;
        for (const cpp_int &q : primes) {
            if (q == 0) continue;
            cpp_int exp = phi / q;
            if (mod_pow(g, exp, p) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return 0;
}

// Генерация случайного cpp_int в диапазоне [1, max-1]
cpp_int random_in_range(const cpp_int &max) {
    // используем mt19937_64 для генерации случайных бит
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    // получим число бит в max
    std::string s = max.convert_to<std::string>();
    // простейшая стратегия: генерируем случайную 64-битную последовательность и собираем
    // Пока это не замечательный крипто-генератор, но для демонстрации хватит.
    cpp_int r = 0;
    uint64_t chunks = 8; // 8*64 = 512 бит (достаточно для больших учебных p)
    for (uint64_t i = 0; i < chunks; ++i) {
        uint64_t v = gen();
        r <<= 64;
        r += v;
    }
    r %= (max - 1);
    r += 1;
    return r;
}

// BD-протокол: вход p,g,n — возвращает общий ключ и печатает шаги.
cpp_int burmester_desmedt(const cpp_int &p, const cpp_int &g, int n, bool show_steps = true) {
    // 1) каждый выбирает секрет x_i
    vector<cpp_int> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = random_in_range(p - 1); // секрет в [1, p-2]
    }

    // 2) вычисляем публичные z_i = g^{x_i} mod p и рассылаем
    vector<cpp_int> z(n);
    for (int i = 0; i < n; ++i) {
        z[i] = mod_pow(g, x[i], p);
    }

    if (show_steps) {
        cout << "Параметры:\n p = " << p << "\n g = " << g << "\n n = " << n << "\n\n";
        for (int i = 0; i < n; ++i) {
            cout << "U" << (i+1) << ": x_" << (i+1) << " = " << x[i] << "  z_" << (i+1) << " = " << z[i] << "\n";
        }
        cout << "\n";
    }

    // В BD часто вычисляют s_{i,i+1} = z_{i+1}^{x_i} = g^{x_i x_{i+1}}
    vector<cpp_int> s(n);
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1) % n;
        s[i] = mod_pow(z[nxt], x[i], p);
        if (show_steps) {
            cout << "s_" << (i+1) << ("," ) << (nxt+1) << " = g^{x_" << (i+1) << " x_" << (nxt+1) << "} = " << s[i] << "\n";
        }
    }

    // общий ключ K = prod_{i=1..n} s_{i,i+1} mod p = g^{sum x_i x_{i+1}}
    cpp_int K = 1;
    for (int i = 0; i < n; ++i) {
        K = (K * s[i]) % p;
    }

    if (show_steps) {
        cout << "\nОбщий ключ (из произведения s_{i,i+1}): K = " << K << "\n";
    }

    // проверим альтернативной формулой: возвести g^{sum(x_i x_{i+1})}
    cpp_int sum = 0;
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1) % n;
        sum += x[i] * x[nxt];
    }
    cpp_int K2 = mod_pow(g, sum, p);
    if (show_steps) {
        cout << "Проверка: sum = " << sum << "\n";
        cout << "g^{sum} mod p = " << K2 << "\n";
        if (K == K2) cout << "Проверка пройдена: K == g^{sum}\n";
        else cout << "Ошибка: K != g^{sum}\n";
    }
    return K;
}

int main() {
    // Два режима: демонстрационный для маленьких p (по умолчанию),
    // и пример с вводом своих p,g (если хочется).
    cout << "Burmester-Desmedt demo (C++).\n\n";

    // Демонстрация на малом примере (как в диалоге): p=23, g=5, n=3
    {
        cpp_int p = 23;
        cpp_int g = 5;
        cout << "=== Демонстрация (p=23, g=5, n=3) ===\n";
        burmester_desmedt(p, g, 3, true);
        cout << "====================================\n\n";
    }

    // Пример: найдем генератор для входного небольшого простого p
    {
        cpp_int p;
        cout << "Введите простое p для поиска генератора g (или 0 чтобы пропустить): ";
        std::string sp;
        std::cin >> sp;
        try {
            p = 0;
            for (char c : sp) if (isdigit(c)) { p *= 10; p += (c - '0'); }
        } catch(...) { p = 0; }
        if (p != 0) {
            cout << "Ищем генератор для p = " << p << " (факторизуем p-1 простым пробным делением)...\n";
            cpp_int g = find_generator(p);
            if (g == 0) {
                cout << "Генератор не найден (возможно p слишком велико для пробного деления или p не простое).\n";
            } else {
                cout << "Найден g = " << g << "\n";
                int n;
                cout << "Введите количество участников n (>=3): ";
                std::cin >> n;
                if (n < 2) n = 3;
                burmester_desmedt(p, g, n, true);
            }
        } else {
            cout << "Пропуск ввода p.\n";
        }
    }

    cout << "Готово.\n";
    return 0;
}
