#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdint>
#include <chrono>
using namespace std;

struct point {
    int64_t x;
    int64_t y;
    bool is_infinity;
    point() : x(0), y(0), is_infinity(true) {}
    point(int64_t x, int64_t y) : x(x), y(y), is_infinity(false) {}
};

struct point_hash {
    size_t operator()(const point& p) const noexcept {
        if (p.is_infinity) return std::hash<int64_t>()(-1);
        
        size_t h1 = std::hash<int64_t>()(p.x);
        size_t h2 = std::hash<int64_t>()(p.y);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
    }
};

bool operator==(const point& a, const point& b) noexcept {
    if (a.is_infinity && b.is_infinity) return true;
    if (a.is_infinity != b.is_infinity) return false;
    return a.x == b.x && a.y == b.y;
}

class ECDH {
private:
    int64_t A; 
    int64_t B; 
    int64_t P; 
    std::mt19937_64 rng;

public:
    ECDH(int64_t a, int64_t b, int64_t p) : A(a), B(b), P(p) {
        std::random_device rd;
        rng.seed( (uint64_t(rd()) << 32) ^ uint64_t(rd()) ^ uint64_t(std::chrono::high_resolution_clock::now().time_since_epoch().count()) );
    }

    
    int64_t mod(int64_t x) const {
        int64_t r = x % P;
        if (r < 0) r += P;
        return r;
    }

    
    int64_t mod_pow(int64_t base, int64_t exp) const {
        base = mod(base);
        int64_t res = 1;
        while (exp > 0) {
            if (exp & 1) res = (__int128)res * base % P;
            base = (__int128)base * base % P;
            exp >>= 1;
        }
        return res;
    }

   
    int64_t modInverse(int64_t x) const {
        x = mod(x);
        if (x == 0) throw runtime_error("modInverse: division by zero");
        return mod_pow(x, P - 2);
    }

    
    bool is_on_curve(const point& pt) const {
        if (pt.is_infinity) return true;
        int64_t lhs = mod((__int128)pt.y * pt.y);
        int64_t rhs = mod((__int128)pt.x * pt.x % P * pt.x + A * pt.x + B);
        return lhs == rhs;
    }

   
    point negate(const point& Pnt) const {
        if (Pnt.is_infinity) return Pnt;
        return point(Pnt.x, mod(-Pnt.y));
    }

    
    point add(const point& P1, const point& P2) const {
        if (P1.is_infinity) return P2;
        if (P2.is_infinity) return P1;
        if (P1.x == P2.x && mod(P1.y + P2.y) == 0) return point(); 

        int64_t lambda;
        if (P1.x == P2.x && P1.y == P2.y) {
           
            if (P1.y == 0) return point();
            int64_t num = mod(3 * (__int128)P1.x * P1.x + A);
            int64_t den = modInverse(mod(2 * P1.y));
            lambda = mod((__int128)num * den);
        } else {
            int64_t num = mod(P2.y - P1.y);
            int64_t den = modInverse(mod(P2.x - P1.x));
            lambda = mod((__int128)num * den);
        }

        int64_t xr = mod((__int128)lambda * lambda - P1.x - P2.x);
        int64_t yr = mod((__int128)lambda * (P1.x - xr) - P1.y);
        return point(xr, yr);
    }

    
    point mul(const point& Pnt, int64_t k) const {
        if (k == 0 || Pnt.is_infinity) return point();
        if (k < 0) return mul(negate(Pnt), -k);

        point R; 
        point Q = Pnt;
        int64_t kk = k;
        while (kk > 0) {
            if (kk & 1) R = add(R, Q);
            Q = add(Q, Q);
            kk >>= 1;
        }
        return R;
    }

   
    
    int64_t find_order(const point& Pnt) {
        if (Pnt.is_infinity) return 1; 

       
        int64_t bound = P + 1 + (int64_t)(2 * floor(sqrt((long double)P)));
        int64_t m = (int64_t)ceil(sqrt((long double)bound));
        
        unordered_map<point,int64_t,point_hash> baby;
        baby.reserve(m+1);
        for (int64_t j = 0; j <= m; ++j) {
            point t = mul(Pnt, j);
            baby.emplace(t, j); 
        }
        
        point mP = mul(Pnt, m);
        point neg_mP = negate(mP);

        point cur; 
        for (int64_t k = 0; k <= m; ++k) {
            auto it = baby.find(cur);
            if (it != baby.end()) {
                int64_t j = it->second;
                int64_t candidate = j + k * m;
                if (candidate > 0) {
                    point test = mul(Pnt, candidate);
                    if (test.is_infinity) return candidate; 
                }
            }
            
            cur = add(cur, neg_mP);
        }
        return -1;
    }

    
    vector<point> list_points(bool print = true) {
        vector<point> pts;
        for (int64_t x = 0; x < P; ++x) {
            int64_t rhs = mod((__int128)x*x%P*x + A*x + B);
            for (int64_t y = 0; y < P; ++y) {
                if (mod(y*y) == rhs) {
                    pts.emplace_back(x,y);
                }
            }
        }
        if (print) {
            for (auto &pt : pts) {
                cout << "P=(" << pt.x << "," << pt.y << ") ord=" << find_order(pt) << "\n";
            }
        }
        return pts;
    }

    
    bool is_generator_of_order(const point& G, int64_t r) {
        if (!is_on_curve(G) || G.is_infinity) return false;
        if (mul(G, r).is_infinity == false) return false;
        
        int64_t ord = find_order(G);
        return ord == r;
    }

    
    bool choose_random_generator_of_order(point &outG, int64_t r) {
        
        vector<point> pts = list_points(false);
        
        vector<point> candidates;
        for (auto &pt : pts) {
            if (pt.is_infinity) continue;
            int64_t ord = find_order(pt);
            if (ord == r) candidates.push_back(pt);
        }
        if (candidates.empty()) return false;
        
        std::uniform_int_distribution<size_t> dist(0, candidates.size()-1);
        outG = candidates[dist(rng)];
        return true;
    }

    
    int64_t generate_private(int64_t r) {
        if (r <= 2) throw runtime_error("order too small");
        std::uniform_int_distribution<int64_t> dist(1, r-1);
        return dist(rng);
    }

    
    static void print_point(const point& Pnt) {
        if (Pnt.is_infinity) cout << "O";
        else cout << "(" << Pnt.x << "," << Pnt.y << ")";
    }
    point multi_party_ecd(ECDH& curve, const point& G, const std::vector<int>& priv_keys) {
    
        point shared = G;
    
        for (int k : priv_keys) {
            shared = curve.mul(shared, k);
        }
        return shared;
    }
};


int main() {
    
    ECDH curve(2,2,17);

   
    cout << "Все точки кривой:\n";
    auto pts = curve.list_points();

    
    point G = pts.empty() ? point() : pts.back();
    int64_t order = curve.find_order(G);
    cout << "\nВыбранная точка: "; ECDH::print_point(G); cout << " порядок=" << order << "\n";

    int n;
    std::cout << "Количество участников: ";
    std::cin >> n;

    std::vector<int> priv_keys(n);
    std::vector<point> pub_keys(n);

    // Генерация ключей
    for (int i = 0; i < n; i++) {
        priv_keys[i] = curve.generate_private(order);
        pub_keys[i] = curve.mul(G, priv_keys[i]);
        std::cout << "Участник " << i+1 << " priv = " << priv_keys[i]
                  << "  pub = (" << pub_keys[i].x << ", " << pub_keys[i].y << ")\n";
    }

    // Вычисляем общий секрет
    point shared_secret = curve.multi_party_ecd(curve, G, priv_keys);
    std::cout << "\nОбщий секрет: (" << shared_secret.x << ", " << shared_secret.y << ")\n";
}