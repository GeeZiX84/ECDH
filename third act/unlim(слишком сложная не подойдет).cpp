
#include <bits/stdc++.h>
#include <nlohmann/json.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using json = nlohmann::json;
using boost::multiprecision::cpp_int;
using namespace std;


cpp_int mod_general(cpp_int x, const cpp_int &m) {
    x %= m;
    if (x < 0) x += m;
    return x;
}

cpp_int mod_add_general(const cpp_int &a, const cpp_int &b, const cpp_int &m) {
    return mod_general(a + b, m);
}
cpp_int mod_sub_general(const cpp_int &a, const cpp_int &b, const cpp_int &m) {
    return mod_general(a - b, m);
}
cpp_int mod_mul_general(const cpp_int &a, const cpp_int &b, const cpp_int &m) {
    return mod_general(a * b, m);
}
cpp_int mod_pow_general(cpp_int base, cpp_int exp, const cpp_int &m) {
    base = mod_general(base, m);
    cpp_int res = 1;
    while (exp > 0) {
        if ((exp & 1) != 0) res = mod_general(res * base, m);
        base = mod_general(base * base, m);
        exp >>= 1;
    }
    return res;
}

// extended gcd inverse modulo m (works for arbitrary cpp_int)
bool egcd_inv(const cpp_int &a_in, const cpp_int &m, cpp_int &out_inv) {
    cpp_int a = a_in, r = m;
    cpp_int old_r = a;
    cpp_int s = 0, old_s = 1;
    while (r != 0) {
        cpp_int q = old_r / r;
        cpp_int tmp = old_r - q * r; old_r = r; r = tmp;
        tmp = old_s - q * s; old_s = s; s = tmp;
    }
    // now old_r = gcd(a,m)
    if (old_r != 1 && old_r != -1) return false;
    cpp_int inv = old_s % m;
    if (inv < 0) inv += m;
    out_inv = inv;
    return true;
}

bool is_prime_small(const cpp_int &n) {
    // используется только для небольших l (мы будем проверять l <= 2^32..64)
    // Попробуем конвертировать в unsigned long long и пробное деление
    try {
        unsigned long long v = n.convert_to<unsigned long long>();
        if (v < 2) return false;
        static const unsigned primes_tab[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
        for (unsigned p: primes_tab) if (v == p) return true;
        for (unsigned p: primes_tab) if (v % p == 0) return false;
        unsigned long long r = floor(sqrt((long double)v));
        for (unsigned long long d = 101; d <= r; d += 2) {
            if (v % d == 0) return false;
        }
        return true;
    } catch(...) {
        // Если l слишком большой для small-primality, для безопасности вернём false (не будем считать Elkies)
        return false;
    }
}

// Лежандров символ (через Эйлера) в поле mod m, где m - простое (ожидается)
// Возвращает: 0 если 0, 1 если квадратный вычет, -1 если невырожденный невчет
int legendre_symbol_mod(const cpp_int &a, const cpp_int &p) {
    cpp_int aa = mod_general(a, p);
    if (aa == 0) return 0;
    cpp_int r = mod_pow_general(aa, (p - 1) / 2, p);
    if (r == 1) return 1;
    if (r == p - 1) return -1;
    // Нечёткий случай (в теории не должен быть)
    return -1;
}

// Tonelli-Shanks для корня по модулю p (p простое, odd)
cpp_int tonelli_shanks_mod(cpp_int n, const cpp_int &p) {
    n = mod_general(n, p);
    if (n == 0) return 0;
    if (p == 2) return n;
    if (mod_pow_general(n, (p - 1) / 2, p) != 1) throw runtime_error("no sqrt");
    if (p % 4 == 3) return mod_pow_general(n, (p + 1) / 4, p);

    // p-1 = q * 2^s
    cpp_int q = p - 1;
    unsigned long long s = 0;
    while ((q & 1) == 0) { q >>= 1; ++s; }

    // find z non-residue
    cpp_int z = 2;
    while (mod_pow_general(z, (p - 1) / 2, p) != p - 1) ++z;

    cpp_int c = mod_pow_general(z, q, p);
    cpp_int x = mod_pow_general(n, (q + 1) / 2, p);
    cpp_int t = mod_pow_general(n, q, p);
    unsigned long long m = s;

    while (t != 1) {
        unsigned long long i = 1;
        cpp_int tt = mod_pow_general(t, 2, p);
        while (tt != 1) {
            tt = mod_pow_general(tt, 2, p);
            ++i;
            if (i == m) throw runtime_error("Tonelli failed");
        }
        unsigned long long shift = m - i - 1;
        cpp_int b = mod_pow_general(c, cpp_int(1ULL) << shift, p);
        x = mod_general(x * b, p);
        c = mod_general(b * b, p);
        t = mod_general(t * c, p);
        m = i;
    }
    return x;
}

//
// Основной класс и структуры
//

struct point {
    cpp_int x;
    cpp_int y;
    bool is_infinity;
    point() : x(0), y(0), is_infinity(true) {}
    point(const cpp_int &xx, const cpp_int &yy) : x(xx), y(yy), is_infinity(false) {}
};
struct point_hash
{
     size_t operator()(const point &p) const noexcept
     {
          if (p.is_infinity)
               return std::hash<long long>()(-1);

          size_t h1 = std::hash<std::string>()(p.x.convert_to<string>());
          size_t h2 = std::hash<std::string>()(p.y.convert_to<string>());
          return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
     }
};

struct cpp_int_hash
{
    size_t operator()(const cpp_int &x) const noexcept
    {
        std::string s;
        try { s = x.convert_to<string>(); } catch(...) { s = "0"; }
        return std::hash<std::string>()(s);
    }
};

bool operator==(const point &a, const point &b) noexcept
{
     if (a.is_infinity && b.is_infinity) return true;
     if (a.is_infinity != b.is_infinity) return false;
     return a.x == b.x && a.y == b.y;
}

std::ostream& operator<<(std::ostream &os, const point &p) {
    if (p.is_infinity) {
        os << "O";
    } else {
        os << p.x << "," << p.y;
    }
    return os;
}

class ECDH {
public:
    cpp_int A;
    cpp_int B;
    cpp_int P;
    std::mt19937_64 rng;

    ECDH(const cpp_int &a, const cpp_int &b, const cpp_int &p) : A(a), B(b), P(p) {
        A = mod(A);
        B = mod(B);
        rng.seed((uint64_t)std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }

    // методы нормальной работы -- все по модулю P
    cpp_int mod(const cpp_int &x) const {
        return mod_general(x, P);
    }
    cpp_int mod_pow(cpp_int base, cpp_int exp) const {
        return mod_pow_general(base, exp, P);
    }
    cpp_int modInverse(const cpp_int &x) const {
        cpp_int inv;
        if (!egcd_inv(mod_general(x, P), P, inv)) throw runtime_error("modInverse: no inverse");
        return inv;
    }

    bool is_on_curve(const point &pt) const {
        if (pt.is_infinity) return true;
        cpp_int lhs = mod(pt.y * pt.y);
        cpp_int rhs = mod(pt.x * pt.x * pt.x + A * pt.x + B);
        return lhs == rhs;
    }

    point negate(const point &Pnt) const {
        if (Pnt.is_infinity) return Pnt;
        return point(Pnt.x, mod(-Pnt.y));
    }

    point add(const point &P1, const point &P2) const {
        if (P1.is_infinity) return P2;
        if (P2.is_infinity) return P1;
        if (P1.x == P2.x) {
            if (mod(P1.y + P2.y) == 0) return point(); // infinity
            // doubling
            if (P1.y == 0) return point();
            cpp_int num = mod(3 * P1.x * P1.x + A);
            cpp_int den = mod(2 * P1.y);
            cpp_int inv = modInverse(den);
            cpp_int lambda = mod(num * inv);
            cpp_int xr = mod(lambda * lambda - P1.x - P2.x);
            cpp_int yr = mod(lambda * (P1.x - xr) - P1.y);
            return point(xr, yr);
        } else {
            cpp_int num = mod(P2.y - P1.y);
            cpp_int den = mod(P2.x - P1.x);
            cpp_int inv = modInverse(den);
            cpp_int lambda = mod(num * inv);
            cpp_int xr = mod(lambda * lambda - P1.x - P2.x);
            cpp_int yr = mod(lambda * (P1.x - xr) - P1.y);
            return point(xr, yr);
        }
    }

    point mul(const point &Pnt, cpp_int k) const {
        if (k == 0 || Pnt.is_infinity) return point();
        if (k < 0) return mul(negate(Pnt), -k);
        point R; // infinity
        point Q = Pnt;
        cpp_int kk = k;
        while (kk > 0) {
            if ((kk & 1) != 0) R = add(R, Q);
            Q = add(Q, Q);
            kk >>= 1;
        }
        return R;
    }

    // find_order via BSGS (медленно для больших p) - оставлен
    cpp_int find_order(const point &Pnt) const {
        if (Pnt.is_infinity) return 1;
        long double pd = 0;
        try { pd = (long double)P.convert_to<long double>(); } catch(...) { pd = 0; }
        long double sq = floor(sqrt(pd));
        cpp_int bound = P + 1 + cpp_int(2 * floor(sq));
        long double bound_ld = (long double) (P.convert_to<long double>() + 1.0L + 2.0L * sqrt((long double)P.convert_to<long double>()));
        cpp_int m = cpp_int( (long long) ceil(sqrt((long double)bound_ld)) );
        if (m <= 0) m = 1;
        unordered_map<string, cpp_int> baby;
        baby.reserve( (size_t) ( (m.convert_to<unsigned long long>()>0) ? m.convert_to<unsigned long long>() : 1 ) );
        point cur;
        for (cpp_int j = 0; j <= m; ++j) {
            string key = (cur.is_infinity ? string("I") : (cur.x.convert_to<string>() + "|" + cur.y.convert_to<string>()));
            baby.emplace(key, j);
            cur = add(cur, Pnt);
        }
        point mP = mul(Pnt, m);
        point giant;
        for (cpp_int k = 1; k <= m + 1; ++k) {
            giant = add(giant, mP);
            string key = (giant.is_infinity ? string("I") : (giant.x.convert_to<string>() + "|" + giant.y.convert_to<string>()));
            auto it = baby.find(key);
            if (it != baby.end()) {
                cpp_int j = it->second;
                cpp_int candidate = k * m - j;
                if (candidate > 0) {
                    if (mul(Pnt, candidate).is_infinity) return candidate;
                }
            }
        }
        return -1;
    }

    // count_points (простой / медленный) - для небольших p
    cpp_int count_points() const {
        cpp_int total = 1;
        bool use_ull = false;
        unsigned long long P_ull = 0;
        try { P_ull = P.convert_to<unsigned long long>(); use_ull = (P == cpp_int(P_ull)); } catch(...) { use_ull = false; }
        if (use_ull) {
            unordered_map<unsigned long long, int> leg_cache;
            for (unsigned long long x = 0; x < P_ull; ++x) {
                cpp_int xx = cpp_int(x);
                cpp_int rhs = mod(xx * xx * xx + A * xx + B);
                if (rhs == 0) { total += 1; continue; }
                int ls = legendre_symbol(rhs);
                if (ls == 1) total += 2;
            }
        } else {
            unordered_map<string, int> leg_cache;
            for (cpp_int x = 0; x < P; ++x) {
                cpp_int rhs = mod(x * x * x + A * x + B);
                if (rhs == 0) { total += 1; continue; }
                int ls = legendre_symbol(rhs);
                if (ls == 1) total += 2;
            }
        }
        return total;
    }

    cpp_int generate_private(const cpp_int &order) {
        static std::mt19937_64 rng_local((uint64_t)std::chrono::high_resolution_clock::now().time_since_epoch().count());
        cpp_int acc = 0;
        for (int i = 0; i < 8; ++i) {
            uint64_t r = rng_local();
            acc <<= 64;
            acc += cpp_int(r);
        }
        if (order > 2) acc = acc % (order - 1) + 1;
        else acc = 1;
        return acc;
    }

    int legendre_symbol(const cpp_int &a) const {
        cpp_int aa = mod_general(a, P);
        if (aa == 0) return 0;
        cpp_int r = mod_pow_general(aa, (P - 1) / 2, P);
        if (r == 1) return 1;
        if (r == P - 1) return -1;
        return -1;
    }

    cpp_int tonelli_shanks(const cpp_int &n) const {
        return tonelli_shanks_mod(n, P);
    }

    point random_point() {
        while (true) {
            cpp_int x = 0;
            for (int i = 0; i < 4; ++i) {
                uint64_t r = rng();
                x <<= 64;
                x += cpp_int(r);
            }
            x = mod(x);
            cpp_int rhs = mod(x * x * x + A * x + B);
            if (legendre_symbol(rhs) != 1) continue;
            cpp_int y = tonelli_shanks(rhs);
            return point(x, y);
        }
    }

    // Проверяет, является ли l простым и является ли он Elkies для кривой.
    // Возвращает true если Elkies, и записывает t mod l в t_mod_l_out.
    bool is_elkies_and_get_t(cpp_int l, cpp_int &t_mod_l_out) const {
        // требуем, чтобы l было небольшим простым (без дорогостоящего MR)
        if (!is_prime_small(l)) return false;
        // считаем #E(F_l) перебором x в [0..l-1]  (работает только для малых l)
        cpp_int points = 1; // бесконечность
        for (cpp_int x = 0; x < l; ++x) {
            cpp_int rhs = mod_general(x * x * x + A * x + B, l);
            if (rhs == 0) { points += 1; continue; }
            cpp_int r = mod_pow_general(rhs, (l - 1) / 2, l);
            if (r == 1) points += 2;
        }
        cpp_int tmod = mod_general((l + 1) - points, l);
        t_mod_l_out = tmod;
        // delta = t^2 - 4p mod l
        cpp_int delta = mod_general(mod_pow_general(tmod, 2, l) - mod_mul_general(4, mod_general(P, l), l), l);
        int leg = legendre_symbol_mod(delta, l);
        return (leg == 1);
    }

}; // class ECDH

//
// CRT
//
cpp_int crt_combine(const vector<cpp_int> &mods, const vector<cpp_int> &rems) {
    if (mods.empty() || mods.size() != rems.size()) throw runtime_error("crt_combine: bad args");
    cpp_int M = 1;
    for (auto &m : mods) M *= m;
    cpp_int x = 0;
    for (size_t i = 0; i < mods.size(); ++i) {
        cpp_int Mi = M / mods[i];
        cpp_int inv;
        if (!egcd_inv(Mi % mods[i], mods[i], inv)) throw runtime_error("crt: no inverse");
        x += rems[i] * Mi * inv;
    }
    x = mod_general(x, M);
    return x;
}

//
// SEA: простой цикл по малым простым l используя только Elkies primes (через is_elkies_and_get_t).
// Собираем t mod (product of l) с CRT.
//
cpp_int sea_compute_trace(const ECDH &E, const vector<int> &small_primes) {
    vector<cpp_int> mods;
    vector<cpp_int> rems;
    cpp_int prod = 1;
    long double bound = 4.0L * sqrt((long double) E.P.convert_to<long double>()); // условие: prod > 4*sqrt(p)
    for (int l_int : small_primes) {
        cpp_int l = cpp_int(l_int);
        if (l == E.P) continue;
        cpp_int t_l;
        bool elk = E.is_elkies_and_get_t(l, t_l);
        if (!elk) continue; // пропускаем Atkin (упрощение)
        mods.push_back(l);
        rems.push_back(t_l);
        prod *= l;
        if ((long double)prod.convert_to<long double>() > bound) break;
    }
    if (mods.empty()) throw runtime_error("SEA: no Elkies residues found");
    cpp_int t = crt_combine(mods, rems);
    // приведение t к диапазону [-2sqrt(p), 2sqrt(p)] возможно позже
    return t;
}

//
// parse helper
//
cpp_int parse_number_to_cppint(const std::string &s_raw) {
    std::string s = s_raw;
    if (s.empty()) throw runtime_error("empty number string");
    bool negative = false;
    size_t pos = 0;
    if (s[0] == '+') pos = 1;
    if (s[pos] == '-') { negative = true; ++pos; }
    std::string body = s.substr(pos);
    cpp_int value = 0;
    if (body.size() >= 2 && (body[0] == '0' && (body[1]=='x' || body[1]=='X'))) {
        for (size_t i = 2; i < body.size(); ++i) {
            char c = body[i];
            int digit = 0;
            if (c >= '0' && c <= '9') digit = c - '0';
            else if (c >= 'a' && c <= 'f') digit = 10 + (c - 'a');
            else if (c >= 'A' && c <= 'F') digit = 10 + (c - 'A');
            else throw runtime_error("invalid hex char in " + body);
            value <<= 4;
            value += digit;
        }
    } else {
        for (char c : body) {
            if (c < '0' || c > '9') throw runtime_error("invalid decimal char in " + body);
            value *= 10;
            value += (c - '0');
        }
    }
    if (negative) value = -value;
    return value;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream f("test2.json");
    if (!f.is_open()) { cerr << "Не удалось открыть test2.json\n"; return 1; }
    json data; f >> data;

    // Набор небольших простых для SEA (можно расширять)
    vector<int> small_primes = {2,3,5,7,11,13,17,19,23,29,31,37,41};

    for (size_t i = 0; i < data.size(); ++i) {
        try {
            cpp_int p = parse_number_to_cppint(data[i]["field"]["p"].get<string>());
            cpp_int a = parse_number_to_cppint(data[i]["a"].get<string>());
            cpp_int b = parse_number_to_cppint(data[i]["b"].get<string>());

            ECDH ecdh(a, b, p);
            cout << "Curve #" << (i+1) << " p=" << p << " a=" << a << " b=" << b << "\n";

            // Если p небольшое — можно просто сосчитать точки (быстро)
            bool p_small = false;
            try { unsigned long long pu = p.convert_to<unsigned long long>(); p_small = true; } catch(...) { p_small = false; }

            cpp_int group_order;
            if (p_small && p < 1000000) {
                cout << "  small p -> brute count\n";
                group_order = ecdh.count_points();
                cout << "  group order = " << group_order << "\n";
            } else {
                cout << "  Running SEA (Elkies-only simplification). This may still be heavy for large p.\n";
                try {
                    cpp_int t_mod_prod = sea_compute_trace(ecdh, small_primes);
                    // восстановление t из t_mod_prod не полностью корректно без нормализации, но получим предположительный t:
                    // Мы получили t modulo M, где M = product(small_primes_used).
                    // В простом варианте ограничимся тем, что выведем p+1 - t_mod_prod (приближение).
                    cpp_int t = t_mod_prod; // здесь нужен CRT/нормализация в диапазоне [-2sqrt(p),2sqrt(p)]
                    group_order = p + 1 - t;
                    cout << "  (approx) trace t (mod product) = " << t << "\n";
                    cout << "  approx group order = " << group_order << "\n";
                } catch (const exception &ex) {
                    cerr << "  SEA failed/insufficient: " << ex.what() << "\n";
                    cout << "  Falling back to brute count (if possible)\n";
                    group_order = ecdh.count_points();
                    cout << "  group order = " << group_order << "\n";
                }
            }

            // Поиск генератора (искать точку порядка == group_order)
            point G;
            bool found = false;
            const int MAX_TRIALS = 2000;
            for (int ttry = 0; ttry < MAX_TRIALS && !found; ++ttry) {
                point cand = ecdh.random_point();
                // проверяем ord(cand) == group_order (дорого). Лучше проверять cand * group_order == O
                point Q = ecdh.mul(cand, group_order);
                if (Q.is_infinity) {
                    G = cand;
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "  Warning: generator not found after " << MAX_TRIALS << " trials.\n";
                continue;
            }
            cout << "  Found G = " << G << "\n";

            cpp_int privA = ecdh.generate_private(group_order);
            cpp_int privB = ecdh.generate_private(group_order);
            point pubA = ecdh.mul(G, privA);
            point pubB = ecdh.mul(G, privB);
            point secA = ecdh.mul(pubB, privA);
            point secB = ecdh.mul(pubA, privB);

            if (secA == secB && !secA.is_infinity) {
                cout << "  Shared secret = " << secA << "\n";
                cout << "Success: shared secrets match\n";
            } else {
                cout << "Error: shared secrets do not match or are infinity\n";
            }

        } catch (const exception &ex) {
            cerr << "Ошибка при обработке кривой #" << (i+1) << ": " << ex.what() << "\n";
        }
    }

    return 0;
}
