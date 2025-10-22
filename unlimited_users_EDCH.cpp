#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <functional>
#include <nlohmann/json.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>

using boost::multiprecision::cpp_int;
using json = nlohmann::json;
using namespace std;
struct point
{
     cpp_int x;
     cpp_int y;
     bool is_infinity;
     point() : x(0), y(0), is_infinity(true) {}
     point(cpp_int x, cpp_int y) : x(x), y(y), is_infinity(false) {}
     size_t size() const noexcept {
         return sizeof(point);
     }
};


struct point_hash
{
     size_t operator()(const point &p) const noexcept
     {
          if (p.is_infinity)
               return std::hash<cpp_int>()(-1);

          size_t h1 = std::hash<cpp_int>()(p.x);
          size_t h2 = std::hash<cpp_int>()(p.y);
          return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
     }
};

struct int64_hash
{
    size_t operator()(const cpp_int &x) const noexcept
    {
        // Базовый хэш от числа
        size_t h1 = std::hash<cpp_int>()(x);

        // Немного "перемешаем" для лучшего распределения
        return h1 ^ (0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

std::ostream& operator<<(std::ostream& os, const point& p)
{
    if (p.is_infinity)
    {
        os << " ";
    }
    else
    {
        os << p.x << "," << p.y << ",";
    }
    return os;
}


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
     cpp_int A;
     cpp_int B;
     cpp_int P;
     std::mt19937_64 rng;

public:
     ECDH(cpp_int a, cpp_int b, cpp_int p) : A(a), B(b), P(p){}
     
     

     
     cpp_int mod(cpp_int x) const
     {
          cpp_int order = x % P;
          if (order < 0)
               order += P;
          return order;
     }

     cpp_int mod_pow(cpp_int base, cpp_int exp) const
     {
          base = mod(base);
          cpp_int res = 1;
          while (exp > 0)
          {
               if (exp & 1)
                    res = (cpp_int)res * base % P;
               base = (cpp_int)base * base % P;
               exp >>= 1;
          }
          return res;
     }

     cpp_int modInverse(cpp_int x) const
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
          cpp_int lhs = mod((cpp_int)pt.y * pt.y);
          cpp_int rhs = mod((cpp_int)pt.x * pt.x % P * pt.x + A * pt.x + B);
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

          cpp_int lambda;
          if (P1.x == P2.x && P1.y == P2.y)
          {

               if (P1.y == 0)
                    return point();
               cpp_int num = mod(3 * (cpp_int)P1.x * P1.x + A);
               cpp_int den = modInverse(mod(2 * P1.y));
               lambda = mod((cpp_int)num * den);
          }
          else
          {
               cpp_int num = mod(P2.y - P1.y);
               cpp_int den = modInverse(mod(P2.x - P1.x));
               lambda = mod((cpp_int)num * den);
          }

          cpp_int xr = mod((cpp_int)lambda * lambda - P1.x - P2.x);
          cpp_int yr = mod((cpp_int)lambda * (P1.x - xr) - P1.y);
          return point(xr, yr);
     }

     point mul(const point &Pnt, cpp_int k) const
     {
          if (k == 0 || Pnt.is_infinity)
               return point();
          if (k < 0)
               return mul(negate(Pnt), -k);

          point R;
          point Q = Pnt;
          cpp_int kk = k;
          while (kk > 0)
          {
               if (kk & 1)
                    R = add(R, Q);
               Q = add(Q, Q);
               kk >>= 1;
          }
          return R;
     }

     cpp_int find_order(const point &Pnt)
     {
          if (Pnt.is_infinity)
               return 1;

          cpp_int bound = P + 1 + (cpp_int)(2 * floor(sqrt((long double)P)));
          cpp_int m = (cpp_int)ceil(sqrt((long double)bound));

          unordered_map<point, cpp_int, point_hash> baby;
          baby.reserve(int64_t(m) + 1);

          // baby steps
          point cur;
          for (cpp_int j = 0; j <= m; ++j)
          {
               baby.emplace(cur, j);
               cur = add(cur, Pnt);
          }

          // giant steps
          point mP = mul(Pnt, m);
          point giant;
          for (cpp_int k = 1; k <= m + 1; ++k)
          {
               giant = add(giant, mP);
               auto it = baby.find(giant);
               if (it != baby.end())
               {
                    cpp_int j = it->second;
                    cpp_int candidate = k * m - j;
                    if (candidate > 0)
                    {
                         if (mul(Pnt, candidate).is_infinity)
                              return candidate;
                    }
               }
          }
          return -1;
     }
     

     vector<point> list_points(bool print = false)
     {
          vector<point> pts;

          for (cpp_int x = 0; x < P; ++x)
          {
               cpp_int rhs = mod(((x * x) % P * x + A * x + B) % P);

               // Проверяем, является ли rhs квадратичным вычетом по модулю P
               if (rhs == 0)
               {
                    pts.emplace_back(x, 0);
                    continue;
               }

               if (mod_pow(rhs, (P - 1) / 2) == 1) // символ Лежандра == 1
               {
                    // На малом P можно просто найти y перебором
                    for (cpp_int y = 1; y < P; ++y)
                    {
                         if (mod(y * y) == rhs)
                         {
                              pts.emplace_back(x, y);
                              if (y != 0)
                              pts.emplace_back(x, P - y); // второй корень
                              break;
                         }
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


     bool is_generator_of_order(const point &G, cpp_int order)
     {
          if (!is_on_curve(G) || G.is_infinity)
               return false;
          if (mul(G, order).is_infinity == false)
               return false;

          cpp_int ord = find_order(G);
          return ord == order;
     }

     bool choose_random_generator_of_order(point &outG, cpp_int order)
     {

          vector<point> pts = list_points(false);

          vector<point> candidates;
          for (auto &pt : pts)
          {
               if (pt.is_infinity)
                    continue;
               cpp_int ord = find_order(pt);
               if (ord == order)
                    candidates.push_back(pt);
          }
          if (candidates.empty())
               return false;

          std::uniform_int_distribution<size_t> dist(0, candidates.size() - 1);
          outG = candidates[dist(rng)];
          return true;
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

     static void print_point(const point &Pnt)
     {
          if (Pnt.is_infinity)
               cout << "O";
          else
               cout << "(" << Pnt.x << "," << Pnt.y << ")";
     }
};
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
// 1) Ring protocol (последовательный обмен).
//    В локальной имитации: для каждого i стартуем с pub[i] = a_i * G и многократно возводим в приват соседей.
//    В реальной сети: участник i получает X от (i-1), вычисляет X^{a_i} и передаёт дальше.
//    Результат: для любого i в конце получится S = (prod a_j) * G.
point multi_ring_shared(const ECDH &curve, const point &G, const vector<cpp_int> &priv, const vector<point> &pub) {
    int n = (int)priv.size();
    if (n == 0) return point();
    // моделируем локально: для i=0..n-1 считаем секрет, проверим совпадение
    // Вернём просто секрет для i=0 (все должны совпадать)
    point cur = pub[0]; // g^{a0}
    for (int step = 1; step <= n-1; ++step) {
        int idx = (0 + step) % n;
        // cur = cur^{a_idx}  (в эллиптической формулировке — умножение точки на скаляр)
        cur = curve.mul(cur, priv[idx]);
    }
    return cur;
}

// 2) Direct local: product of privs mod order, then mul(G, product)
//    Требует: знание порядка группы (order > 0)
point multi_product_shared(const ECDH &curve, const point &G, const vector<cpp_int> &priv, const cpp_int &order) {
    if (order <= 0) throw runtime_error("order <= 0 required for product method");
    cpp_int prod = 1;
    for (auto &a : priv) {
        prod = (prod * (a % order)) % order;
    }
    // итоговый секрет = prod * G
    return curve.mul(G, prod);
}

int main()
{
     std::ofstream out;
     out.open("D:/GitHub/ECDH/curves4.txt", ios::app); 

     out << "A," << "B," << "P,"
     << "time_gen(micros),time_prod(micros),"<<"shared_prod_x,shared_prod_y,"<< "shared_ring_x,shared_ring_y," << "order," << " time(ms)" << endl;
     out.close();
     ifstream f("testbd.json");
     if (!f.is_open()) { cerr << "Не удалось открыть test2.json\n"; return 1; }
     json data; f >> data;
     auto start2 = chrono::high_resolution_clock::now();
     for (size_t i = 0; i < data.size(); ++i) {
          cpp_int p = parse_number_to_cppint(data[i]["field"]["p"].get<string>());
          cpp_int a = parse_number_to_cppint(data[i]["a"].get<string>());
          cpp_int b = parse_number_to_cppint(data[i]["b"].get<string>());    
          auto start = chrono::high_resolution_clock::now(); 
          ECDH curve(a, b, p);
          point G;
          cpp_int order = -1;
          auto pts = curve.list_points(false);
          for (auto &pt : pts) {
               if (pt.is_infinity) continue;
               cpp_int ord = curve.find_order(pt);
               if (ord > 1 and curve.is_generator_of_order(pt, ord) == true) {
                    G = pt;
                    order = ord;
                    break;
               }
          }
          if (order <= 1) {
               cerr << "Не удалось найти точку-генератор для этой кривой (или P слишком велик для brute list).\n";
               continue;
          }
          cout << "Найдена точка G=" << G << " order=" << order << "\n";

          // читаем количество пользователей
          cout << "Введите количество пользователей (положительное целое): ";
          int n = 4;
          // Засекаем время генерации приватных и публичных ключей
          vector<cpp_int> priv(n);
          vector<point> pub(n);

          auto t_gen_start = chrono::high_resolution_clock::now();
          for (int i = 0; i < n; ++i) {
               priv[i] = curve.generate_private(order); // используем order если известен
               pub[i] = curve.mul(G, priv[i]);
          }

          point shared_ring = multi_ring_shared(curve, G, priv, pub);
          auto t_gen_end = chrono::high_resolution_clock::now();
          auto t_gen_start2 = chrono::high_resolution_clock::now();
          bool have_order = (order > 0);
          point shared_prod;
          
          if (have_order) {
               
               shared_prod = multi_product_shared(curve, G, priv, order);
               
               
          }
          auto t_gen_end2 = chrono::high_resolution_clock::now();
          auto time_gen = chrono::duration_cast<chrono::microseconds>(t_gen_end - t_gen_start).count();
          cout << "Общий секрет (умножением G на сумму приватных): " << shared_prod << "\n";
          cout << "Время генерации ключей (micros): " << time_gen << "\n";
          auto time_prod = chrono::duration_cast<chrono::microseconds>(t_gen_end2 - t_gen_start2).count();
          cout << "Время вычисления общего секрета product (micros): " << chrono::duration_cast<chrono::microseconds>(t_gen_end2 - t_gen_start2).count() << "\n";
          cout << "Проверка совпадения секретов: ";
          if (shared_ring == shared_prod){ 
               cout << "совпадают\n";
          }
          else{
               cout << "НЕ совпадают\n";
               continue;
          }
          cout << "Общий секрет (последовательным возведением в приват соседей): " << shared_ring << "\n";
                       
          
          std::ofstream out;          // поток для записи
          out.open("D:/GitHub/ECDH/curves4.txt", ios::app); 
               if (out.is_open()){
                    auto stop = chrono::high_resolution_clock::now();
                    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
                                  
                    out << a << "," << b << "," << p << "," << time_gen << "," << time_prod << ","
                    <<shared_prod << shared_ring.size() << "," << order << "," << duration.count() << endl;
                    out.close();
               }
     }              
     auto stop2 = chrono::high_resolution_clock::now();
     auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
     
     cout << "Total time taken: " << duration2.count() << " milliseconds" << endl;
     
     return 0;
}