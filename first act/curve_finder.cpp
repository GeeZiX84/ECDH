#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdint>
#include <chrono>
using namespace std;
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>


struct point
{
     int64_t x;
     int64_t y;
     bool is_infinity;
     point() : x(0), y(0), is_infinity(true) {}
     point(int64_t x, int64_t y) : x(x), y(y), is_infinity(false) {}
};

struct point_hash
{
     size_t operator()(const point &p) const noexcept
     {
          if (p.is_infinity)
               return std::hash<int64_t>()(-1);

          size_t h1 = std::hash<int64_t>()(p.x);
          size_t h2 = std::hash<int64_t>()(p.y);
          return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
     }
};

bool operator==(const point &a, const point &b) noexcept
{
     if (a.is_infinity && b.is_infinity)
          return true;
     if (a.is_infinity != b.is_infinity)
          return false;
     return a.x == b.x && a.y == b.y;
}

class ECDH
{
private:
     int64_t A;
     int64_t B;
     int64_t P;
     std::mt19937_64 rng;

public:
     ECDH(int64_t a, int64_t b, int64_t p) : A(a), B(b), P(p){}
     
     int64_t mod(int64_t x) const
     {
          int64_t order = x % P;
          if (order < 0)
               order += P;
          return order;
     }

     int64_t mod_pow(int64_t base, int64_t exp) const
     {
          base = mod(base);
          int64_t res = 1;
          while (exp > 0)
          {
               if (exp & 1)
                    res = (int64_t)res * base % P;
               base = (int64_t)base * base % P;
               exp >>= 1;
          }
          return res;
     }

     int64_t modInverse(int64_t x) const
     {
          x = mod(x);
          if (x == 0)
               throw runtime_error("modInverse: division by zero");
          return mod_pow(x, P - 2);
     }

     bool is_on_curve(const point &pt) const
     {
          if (pt.is_infinity)
               return true;
          int64_t lhs = mod((int64_t)pt.y * pt.y);
          int64_t rhs = mod((int64_t)pt.x * pt.x % P * pt.x + A * pt.x + B);
          return lhs == rhs;
     }

     point negate(const point &Pnt) const
     {
          if (Pnt.is_infinity)
               return Pnt;
          return point(Pnt.x, mod(-Pnt.y));
     }

     point add(const point &P1, const point &P2) const
     {
          if (P1.is_infinity)
               return P2;
          if (P2.is_infinity)
               return P1;
          if (P1.x == P2.x && mod(P1.y + P2.y) == 0)
               return point();

          int64_t lambda;
          if (P1.x == P2.x && P1.y == P2.y)
          {

               if (P1.y == 0)
                    return point();
               int64_t num = mod(3 * (int64_t)P1.x * P1.x + A);
               int64_t den = modInverse(mod(2 * P1.y));
               lambda = mod((int64_t)num * den);
          }
          else
          {
               int64_t num = mod(P2.y - P1.y);
               int64_t den = modInverse(mod(P2.x - P1.x));
               lambda = mod((int64_t)num * den);
          }

          int64_t xr = mod((int64_t)lambda * lambda - P1.x - P2.x);
          int64_t yr = mod((int64_t)lambda * (P1.x - xr) - P1.y);
          return point(xr, yr);
     }

     point mul(const point &Pnt, int64_t k) const
     {
          if (k == 0 || Pnt.is_infinity)
               return point();
          if (k < 0)
               return mul(negate(Pnt), -k);

          point R;
          point Q = Pnt;
          int64_t kk = k;
          while (kk > 0)
          {
               if (kk & 1)
                    R = add(R, Q);
               Q = add(Q, Q);
               kk >>= 1;
          }
          return R;
     }

     int64_t find_order(const point &Pnt)
     {
          if (Pnt.is_infinity)
               return 1;

          int64_t bound = P + 1 + (int64_t)(2 * floor(sqrt((long double)P)));
          int64_t m = (int64_t)ceil(sqrt((long double)bound));
          
          unordered_map<point, int64_t, point_hash> baby;
          baby.reserve(m + 1);
          for (int64_t j = 0; j <= m; ++j)
          {
               point t = mul(Pnt, j);
               baby.emplace(t, j);
          }

          point mP = mul(Pnt, m);
          point neg_mP = negate(mP);

          point cur;
          for (int64_t k = 0; k <= m; ++k)
          {
               auto it = baby.find(cur);
               if (it != baby.end())
               {
                    int64_t j = it->second;
                    int64_t candidate = j + k * m;
                    if (candidate > 0)
                    {
                         point test = mul(Pnt, candidate);
                         if (test.is_infinity)
                              return candidate;
                    }
               }

               cur = add(cur, neg_mP);
          }
          return -1;
     }

     vector<point> list_points(bool print = false)
     {
          vector<point> pts;
          for (int64_t x = 0; x < P; ++x)
          {
               int64_t rhs = mod((int64_t)x * x % P * x + A * x + B);
               for (int64_t y = 0; y < P; ++y)
               {
                    if (mod(y * y) == rhs)
                    {
                         pts.emplace_back(x, y);
                    }
               }
          }
          if (print)
          {
               for (auto &pt : pts)
               {
                    cout << "P=(" << pt.x << "," << pt.y << ") ord=" << find_order(pt) << "\n";
               }
          }
          return pts;
     }

     bool is_generator_of_order(const point &G, int64_t order)
     {
          if (!is_on_curve(G) || G.is_infinity)
               return false;
          if (mul(G, order).is_infinity == false)
               return false;

          int64_t ord = find_order(G);
          return ord == order;
     }

     bool choose_random_generator_of_order(point &outG, int64_t order)
     {

          vector<point> pts = list_points(false);

          vector<point> candidates;
          for (auto &pt : pts)
          {
               if (pt.is_infinity)
                    continue;
               int64_t ord = find_order(pt);
               if (ord == order)
                    candidates.push_back(pt);
          }
          if (candidates.empty())
               return false;

          std::uniform_int_distribution<size_t> dist(0, candidates.size() - 1);
          outG = candidates[dist(rng)];
          return true;
     }

     int64_t generate_private(int64_t order)
     {
          if (order <= 2)
               throw runtime_error("order too small");
          std::uniform_int_distribution<int64_t> dist(1, order - 1);
          return dist(rng);
     }

     static void print_point(const point &Pnt)
     {
          if (Pnt.is_infinity)
               cout << "O";
          else
               cout << "(" << Pnt.x << "," << Pnt.y << ")";
     }
};

int main()
{
     auto start2 = chrono::high_resolution_clock::now();
     for (int64_t A = 0; A < 100; A++){
          for (int64_t B = 0; B < 100; B++){
               if (4 * A * A * A + 27 * B * B != 0){
                    auto start = chrono::high_resolution_clock::now();
                    int64_t P = 2267; 
                    ECDH curve(A, B, P);
                    
                    auto pts = curve.list_points();

                    point G = pts.empty() ? point() : pts.back();
                    int64_t order = curve.find_order(G);
                    

                    int64_t da = curve.generate_private(order);
                    int64_t db = curve.generate_private(order);

                    

                    // публичные ключи
                    point Qa = curve.mul(G, da);
                    point Qb = curve.mul(G, db);
                   

                    // общий секрет
                    point Sa = curve.mul(Qb, da);
                    point Sb = curve.mul(Qa, db);
                    
                    std::ofstream out;          // поток для записи
                    out.open("curves.txt", ios::app); 
                    if (Sa == Sb){
                         cout << "Shared secret совпадает.\n";
                         if (out.is_open()){
                              auto stop = chrono::high_resolution_clock::now();
                              auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
                               
                              out << "A=" << A << " B=" << B << " P=" << P << " Время=" << duration.count() << endl;
                              out.close();
                         }
                         

                    }
                    else{
                         cout << "Ошибка: секреты не совпадают!\n";
                    }
               
               }
          }
     }
     auto stop2 = chrono::high_resolution_clock::now();
     auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
     cout << "Total time taken: " << duration2.count() << " milliseconds" << endl;
     return 0;
}