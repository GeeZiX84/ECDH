// console app.cpp
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <functional>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using json = nlohmann::json;
using boost::multiprecision::cpp_int;

// ====================== Константы и версия ======================
const string VERSION = "1.0.0";
const string BUILD_DATE = __DATE__ " " __TIME__;

// ====================== Общие утилиты ======================

cpp_int parse_number_to_cppint(const string &s_raw) {
    string s = s_raw;
    if (s.empty()) throw runtime_error("empty number string");
    
    s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
    
    bool negative = false;
    size_t pos = 0;
    if (s[0] == '+') pos = 1;
    if (s[pos] == '-') { negative = true; ++pos; }
    
    string body = s.substr(pos);
    cpp_int value = 0;
    
    if (body.size() >= 2 && (body[0] == '0' && (body[1]=='x' || body[1]=='X'))) {
        for (size_t i = 2; i < body.size(); ++i) {
            char c = body[i];
            int digit = 0;
            if (c >= '0' && c <= '9') digit = c - '0';
            else if (c >= 'a' && c <= 'f') digit = 10 + (c - 'a');
            else if (c >= 'A' && c <= 'F') digit = 10 + (c - 'A');
            else throw runtime_error("invalid hex char in " + body);
            value = (value << 4) + digit;
        }
    } else {
        for (char c : body) {
            if (c < '0' || c > '9') throw runtime_error("invalid decimal char in " + body);
            value = value * 10 + (c - '0');
        }
    }
    return negative ? -value : value;
}

size_t get_key_size(const cpp_int& key) {
    if (key == 0) return 0;
    string bin_str = key.str(2);
    if (bin_str[0] == '-') bin_str = bin_str.substr(1);
    return bin_str.length();
}

cpp_int random_in_range(const cpp_int &max) {
    if (max <= 2) return 1;
    
    static random_device rd;
    static mt19937_64 gen(rd());
    
    size_t bits = get_key_size(max);
    
    if (bits <= 64) {
        uniform_int_distribution<uint64_t> dist(1, (uint64_t)(max - 1));
        return cpp_int(dist(gen));
    }
    
    cpp_int r = 0;
    size_t chunks = (bits + 63) / 64;
    
    for (size_t i = 0; i < chunks; ++i) {
        uint64_t v = gen();
        r = (r << 64) | cpp_int(v);
    }
    
    r %= (max - 1);
    return (r == 0) ? cpp_int(1) : r;
}

cpp_int mod_pow(cpp_int base, cpp_int exp, const cpp_int &mod) {
    if (mod == 1) return 0;
    base %= mod;
    cpp_int result = 1;
    
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

// ====================== Burmester-Desmedt ======================

struct BDStats {
    long long t_x_us = 0, t_z_us = 0, t_s_us = 0, t_K_us = 0;
    long long total_us = 0;
    size_t key_bits = 0;
};

vector<cpp_int> factorize_small(const cpp_int &n) {
    vector<cpp_int> res;
    cpp_int x = n;
    
    if (x % 2 == 0) {
        res.push_back(2);
        while (x % 2 == 0) x /= 2;
    }
    
    uint64_t limit = 1000000;
    for (uint64_t d = 3; d <= limit && x > 1; d += 2) {
        if (x % d == 0) {
            res.push_back(d);
            while (x % d == 0) x /= d;
        }
    }
    
    if (x > 1) res.push_back(x);
    return res;
}

cpp_int find_generator_near_p(const cpp_int &p, int back_steps = 1000) {
    if (p <= 3) return 0;
    
    cpp_int phi = p - 1;
    auto primes = factorize_small(phi);
    cpp_int start = p - 2;
    cpp_int min_g = (p - 2 - back_steps > 2) ? p - 2 - back_steps : cpp_int(2);
    
    for (cpp_int g = start; g >= min_g; --g) {
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

cpp_int burmester_desmedt(const cpp_int &p, const cpp_int &g, int n, 
                         BDStats &stats, bool show_steps = false) {
    using namespace chrono;
    auto t_total_start = high_resolution_clock::now();

    // Генерация секретов
    auto t_x_start = high_resolution_clock::now();
    vector<cpp_int> x(n);
    for (int i = 0; i < n; ++i) x[i] = random_in_range(p - 1);
    stats.t_x_us = duration_cast<microseconds>(
        high_resolution_clock::now() - t_x_start).count();

    // Публичные ключи
    auto t_z_start = high_resolution_clock::now();
    vector<cpp_int> z(n);
    for (int i = 0; i < n; ++i) z[i] = mod_pow(g, x[i], p);
    stats.t_z_us = duration_cast<microseconds>(
        high_resolution_clock::now() - t_z_start).count();

    if (show_steps) {
        cout << "Burmester-Desmedt Protocol\np = " << p << "\ng = " << g 
             << "\nn = " << n << "\n\n";
        for (int i = 0; i < n; ++i)
            cout << "U" << (i+1) << ": x = " << x[i] << "  z = " << z[i] << "\n";
    }

    // Промежуточные значения
    auto t_s_start = high_resolution_clock::now();
    vector<cpp_int> s(n);
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1) % n;
        s[i] = mod_pow(z[nxt], x[i], p);
    }
    stats.t_s_us = duration_cast<microseconds>(
        high_resolution_clock::now() - t_s_start).count();

    // Общий ключ
    auto t_K_start = high_resolution_clock::now();
    cpp_int K = 1;
    for (int i = 0; i < n; ++i) K = (K * s[i]) % p;
    
    cpp_int sum = 0;
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1) % n;
        sum += x[i] * x[nxt];
    }
    cpp_int K2 = mod_pow(g, sum, p);
    stats.t_K_us = duration_cast<microseconds>(
        high_resolution_clock::now() - t_K_start).count();

    stats.total_us = duration_cast<microseconds>(
        high_resolution_clock::now() - t_total_start).count();
    stats.key_bits = get_key_size(K);

    if (show_steps) {
        cout << "\nKey K = " << K 
             << (K == K2 ? " (VERIFIED)\n" : " (ERROR)\n");
    }

    return K;
}

// ====================== ECDH ======================

struct Point { 
    cpp_int x, y; 
    bool is_infinity;
    
    Point() : x(0), y(0), is_infinity(true) {}
    Point(cpp_int x, cpp_int y) : x(x), y(y), is_infinity(false) {}
    
    size_t size() const { 
        if (is_infinity) return 0;
        return get_key_size(x) + get_key_size(y); 
    }
};

struct PointHash {
    size_t operator()(const Point &p) const {
        if (p.is_infinity) return hash<cpp_int>()(-1);
        size_t h1 = hash<cpp_int>()(p.x);
        size_t h2 = hash<cpp_int>()(p.y);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

bool operator==(const Point &a, const Point &b) {
    if (a.is_infinity && b.is_infinity) return true;
    if (a.is_infinity != b.is_infinity) return false;
    return a.x == b.x && a.y == b.y;
}

class EllipticCurve {
    cpp_int A, B, P;
    
public:
    EllipticCurve(cpp_int a, cpp_int b, cpp_int p) : A(a), B(b), P(p) {}
    
    cpp_int mod(cpp_int x) const { 
        cpp_int r = x % P; 
        return r < 0 ? r + P : r; 
    }
    
    Point add(const Point &P1, const Point &P2) const {
        if (P1.is_infinity) return P2;
        if (P2.is_infinity) return P1;
        if (P1.x == P2.x && mod(P1.y + P2.y) == 0) return Point();
        
        cpp_int lambda, xr, yr;
        if (P1.x == P2.x && P1.y == P2.y) {
            if (P1.y == 0) return Point();
            cpp_int num = mod(3 * P1.x * P1.x + A);
            cpp_int den = mod_pow(2 * P1.y, P - 2, P);
            lambda = mod(num * den);
        } else {
            cpp_int num = mod(P2.y - P1.y);
            cpp_int den = mod_pow(P2.x - P1.x, P - 2, P);
            lambda = mod(num * den);
        }
        
        xr = mod(lambda * lambda - P1.x - P2.x);
        yr = mod(lambda * (P1.x - xr) - P1.y);
        return Point(xr, yr);
    }
    
    Point mul(const Point &Pnt, cpp_int k) const {
        if (k == 0 || Pnt.is_infinity) return Point();
        if (k < 0) {
            Point neg = Point(Pnt.x, mod(-Pnt.y));
            return mul(neg, -k);
        }
        
        Point R;
        Point Q = Pnt;
        while (k > 0) {
            if (k & 1) R = add(R, Q);
            Q = add(Q, Q);
            k >>= 1;
        }
        return R;
    }
    
    bool is_on_curve(const Point &pt) const {
        if (pt.is_infinity) return true;
        cpp_int lhs = mod(pt.y * pt.y);
        cpp_int rhs = mod(pt.x * pt.x * pt.x + A * pt.x + B);
        return lhs == rhs;
    }
    
    bool find_generator(Point &outG, cpp_int &outOrder) {
        // Упрощенный поиск генератора
        for (cpp_int x = 1; x < min(P, cpp_int(100)); ++x) {
            cpp_int rhs = mod(x * x * x + A * x + B);
            
            // Проверяем, является ли rhs квадратичным вычетом
            if (mod_pow(rhs, (P - 1) / 2, P) == 1) {
                // Ищем y перебором
                for (cpp_int y = 0; y < P; ++y) {
                    if (mod(y * y) == rhs) {
                        Point G(x, y);
                        if (!is_on_curve(G)) continue;
                        
                        // Проверяем порядок (упрощенно)
                        Point R = G;
                        cpp_int order = 1;
                        while (!R.is_infinity && order < 100) {
                            R = add(R, G);
                            order++;
                        }
                        
                        if (!R.is_infinity && order >= 2) {
                            outG = G;
                            outOrder = order;
                            return true;
                        }
                        break; // Пробуем следующее x
                    }
                }
            }
        }
        return false;
    }
    
    cpp_int generate_private_key(const cpp_int &order) {
        static mt19937_64 gen(chrono::system_clock::now().time_since_epoch().count());
        if (order <= 2) return 1;
        
        cpp_int r;
        do {
            r = random_in_range(order);
        } while (r == 0);
        return r;
    }
};

// Убрали неиспользуемый параметр G
Point multi_ring_shared(const EllipticCurve &curve, const vector<cpp_int> &priv, 
                       const vector<Point> &pub) {
    int n = (int)priv.size();
    if (n == 0) return Point();
    Point cur = pub[0];
    for (int step = 1; step < n; ++step) {
        int idx = step % n;
        cur = curve.mul(cur, priv[idx]);
    }
    return cur;
}

Point multi_product_shared(const EllipticCurve &curve, const Point &G,
                          const vector<cpp_int> &priv, const cpp_int &order) {
    if (order <= 0) throw runtime_error("order <= 0 required for product method");
    
    cpp_int prod = 1;
    for (auto &a : priv) prod = (prod * (a % order)) % order;
    return curve.mul(G, prod);
}

// ====================== Основная программа ======================

void print_help() {
    cout << "Crypto Protocols v" << VERSION << " (" << BUILD_DATE << ")\n";
    cout << "Usage: crypto_protocols --protocol PROTOCOL --input FILE [options]\n\n";
    cout << "Protocols:\n";
    cout << "  bd     Burmester-Desmedt (DH-based group key exchange)\n";
    cout << "  ecdh   Multi-party ECDH (Elliptic Curve)\n\n";
    cout << "Options:\n";
    cout << "  -p, --protocol    Protocol to use (required)\n";
    cout << "  -i, --input       Input JSON file with curve parameters (required)\n";
    cout << "  -u, --users       Number of users (default: 4)\n";
    cout << "  -o, --output      Output CSV file (default: results.csv)\n";
    cout << "  -v, --verbose     Show detailed output\n";
    cout << "  -h, --help        Show this help\n";
    cout << "  --version         Show version information\n\n";
    cout << "Examples:\n";
    cout << "  crypto_protocols -p bd -i curves.json -u 4\n";
    cout << "  crypto_protocols -p ecdh -i curves.json -u 5 -o ecdh_results.csv -v\n";
}

void print_version() {
    cout << "Crypto Protocols v" << VERSION << "\n";
    cout << "Build date: " << BUILD_DATE << "\n";
    cout << "Protocols: Burmester-Desmedt, Multi-party ECDH\n";
    cout << "License: MIT\n";
}

int main(int argc, char* argv[]) {
    string protocol, input_file, output_file = "results.csv";
    int num_users = 4;
    bool verbose = false;
    
    if (argc < 2) { 
        print_help(); 
        return 1; 
    }
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") { 
            print_help(); 
            return 0; 
        }
        else if (arg == "--version") { 
            print_version(); 
            return 0; 
        }
        else if (arg == "--protocol" || arg == "-p") {
            if (++i < argc) {
                protocol = argv[i];
                if (protocol != "bd" && protocol != "ecdh") {
                    cerr << "Error: unknown protocol '" << protocol << "'\n";
                    return 1;
                }
            }
        } else if (arg == "--input" || arg == "-i") {
            if (++i < argc) input_file = argv[i];
        } else if (arg == "--output" || arg == "-o") {
            if (++i < argc) output_file = argv[i];
        } else if (arg == "--users" || arg == "-u") {
            if (++i < argc) {
                try {
                    num_users = stoi(argv[i]);
                    if (num_users < 2) {
                        cerr << "Error: users must be >= 2\n";
                        return 1;
                    }
                } catch (...) {
                    cerr << "Error: invalid number of users\n";
                    return 1;
                }
            }
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = true;
        }
    }
    
    if (protocol.empty() || input_file.empty()) {
        cerr << "Error: protocol and input file are required\n";
        print_help();
        return 1;
    }
    
    if (verbose) {
        cout << "Starting Crypto Protocols v" << VERSION << "\n";
        cout << "Protocol: " << protocol << "\n";
        cout << "Users: " << num_users << "\n";
        cout << "Input: " << input_file << "\n";
        cout << "Output: " << output_file << "\n\n";
    }
    
    // Загрузка данных
    ifstream in_file(input_file);
    if (!in_file) {
        cerr << "Error: cannot open input file '" << input_file << "'\n";
        return 1;
    }
    
    json data;
    try { 
        in_file >> data; 
    } catch (const exception &e) {
        cerr << "Error parsing JSON: " << e.what() << "\n";
        return 1;
    }
    
    // Выполнение протокола
    ofstream out_file(output_file);
    if (!out_file) {
        cerr << "Error: cannot create output file '" << output_file << "'\n";
        return 1;
    }
    
    auto total_start = chrono::high_resolution_clock::now();
    size_t processed = 0;
    
    try {
        if (protocol == "bd") {
            out_file << "p,g,time_gen_us,time_pub_us,time_inter_us,time_key_us,"
                     << "key_bits,shared_key,total_us\n";
            
            for (size_t i = 0; i < data.size(); ++i) {
                cpp_int p = parse_number_to_cppint(data[i]["field"]["p"].get<string>());
                cpp_int g = find_generator_near_p(p, 100);
                
                if (g == 0) {
                    if (verbose) cerr << "Warning: no generator for p=" << p << "\n";
                    continue;
                }
                
                BDStats stats;
                cpp_int key = burmester_desmedt(p, g, num_users, stats, verbose);
                
                out_file << p << "," << g << "," 
                         << stats.t_x_us << "," << stats.t_z_us << ","
                         << stats.t_s_us << "," << stats.t_K_us << ","
                         << stats.key_bits << "," << key << ","
                         << stats.total_us << "\n";
                
                out_file.flush();
                processed++;
                if (verbose && processed % 10 == 0) {
                    cout << "Processed " << processed << " curves...\n";
                }
            }
            
        } else if (protocol == "ecdh") {
            out_file << "a,b,p,time_ring_us,time_prod_us,key_x,key_y,"
                     << "key_bits,order,total_ms\n";
            
            for (size_t i = 0; i < data.size(); ++i) {
                cpp_int p, a, b;
                try {
                    p = parse_number_to_cppint(data[i]["field"]["p"].get<string>());
                    a = parse_number_to_cppint(data[i]["a"].get<string>());
                    b = parse_number_to_cppint(data[i]["b"].get<string>());
                } catch (const exception &e) {
                    if (verbose) cerr << "Warning: skipping invalid curve data: " << e.what() << "\n";
                    continue;
                }
                
                auto curve_start = chrono::high_resolution_clock::now();
                EllipticCurve curve(a, b, p);
                Point G;
                cpp_int order;
                
                if (!curve.find_generator(G, order)) {
                    if (verbose) cerr << "Warning: no generator for curve " 
                                     << a << "," << b << " mod " << p << "\n";
                    continue;
                }
                
                if (verbose) {
                    cout << "Found generator G=(" << G.x << "," << G.y 
                         << ") order=" << order << "\n";
                }
                
                // Генерация ключей
                vector<cpp_int> priv(num_users);
                vector<Point> pub(num_users);
                
                auto time_start = chrono::high_resolution_clock::now();
                for (int j = 0; j < num_users; ++j) {
                    priv[j] = curve.generate_private_key(order);
                    pub[j] = curve.mul(G, priv[j]);
                }
                auto time_ring = chrono::duration_cast<chrono::microseconds>(
                    chrono::high_resolution_clock::now() - time_start).count();
                
                // Общий секрет
                time_start = chrono::high_resolution_clock::now();
                // Используем исправленную функцию без параметра G
                Point shared_ring = multi_ring_shared(curve, priv, pub);
                Point shared_prod = multi_product_shared(curve, G, priv, order);
                auto time_prod = chrono::duration_cast<chrono::microseconds>(
                    chrono::high_resolution_clock::now() - time_start).count();
                
                // Убрана неиспользуемая переменная total_ms
                auto curve_end = chrono::high_resolution_clock::now();
                auto curve_total_ms = chrono::duration_cast<chrono::milliseconds>(
                    curve_end - curve_start).count();
                
                // Проверка совпадения
                if (!(shared_ring == shared_prod)) {
                    if (verbose) {
                        cerr << "Warning: ring and product results don't match for curve "
                             << a << "," << b << " mod " << p << "\n";
                    }
                }
                
                out_file << a << "," << b << "," << p << ","
                         << time_ring << "," << time_prod << ","
                         << shared_prod.x << "," << shared_prod.y << ","
                         << shared_prod.size() << "," << order << ","
                         << curve_total_ms << "\n";
                
                out_file.flush();
                processed++;
                if (verbose && processed % 10 == 0) {
                    cout << "Processed " << processed << " curves...\n";
                }
            }
        }
        
        out_file.close();
        
        auto total_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::high_resolution_clock::now() - total_start).count();
        
        cout << "\nSUCCESS: Processed " << processed << " curves in " 
             << total_time << " ms\n";
        cout << "Results saved to: " << output_file << "\n";
        
    } catch (const exception &e) {
        cerr << "FATAL ERROR: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}