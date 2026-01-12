#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <iomanip>

#include <openssl/ec.h>
#include <openssl/evp.h>
#include <openssl/bn.h>
#include <openssl/obj_mac.h>
#include <openssl/err.h>

static void ossl_throw(const char* what) {
    unsigned long err = ERR_get_error();
    char buf[256]; ERR_error_string_n(err, buf, sizeof(buf));
    std::string msg = std::string(what) + " | OpenSSL: " + buf;
    throw std::runtime_error(msg);
}

struct OsslCtx {
    BN_CTX* bn = nullptr;
    EC_GROUP* group = nullptr;

    OsslCtx(int nid = NID_X9_62_prime256v1) {
        ERR_load_crypto_strings();
        OpenSSL_add_all_algorithms();
        bn = BN_CTX_new();
        if (!bn) ossl_throw("BN_CTX_new failed");

        group = EC_GROUP_new_by_curve_name(nid);
        if (!group) ossl_throw("EC_GROUP_new_by_curve_name failed");

        EC_GROUP_set_point_conversion_form(group, POINT_CONVERSION_UNCOMPRESSED);
    }
    ~OsslCtx() {
        if (group) EC_GROUP_free(group);
        if (bn) BN_CTX_free(bn);
        EVP_cleanup();
        ERR_free_strings();
    }
};

struct KeyPair {
    BIGNUM* priv = nullptr;       
    EC_POINT* pub = nullptr;       
    KeyPair(): priv(BN_new()), pub(nullptr) {}
    ~KeyPair() { if (priv) BN_free(priv); if (pub) EC_POINT_free(pub); }
};

static KeyPair generate_keypair(const OsslCtx& ctx) {
    KeyPair kp;
    EC_KEY* eckey = EC_KEY_new();
    if (!eckey) ossl_throw("EC_KEY_new failed");
    if (EC_KEY_set_group(eckey, ctx.group) != 1) { EC_KEY_free(eckey); ossl_throw("EC_KEY_set_group failed"); }
    if (EC_KEY_generate_key(eckey) != 1) { EC_KEY_free(eckey); ossl_throw("EC_KEY_generate_key failed"); }

    const BIGNUM* d = EC_KEY_get0_private_key(eckey);
    const EC_POINT* Q = EC_KEY_get0_public_key(eckey);
    if (!d || !Q) { EC_KEY_free(eckey); ossl_throw("EC_KEY_get0_* failed"); }

    if (BN_copy(kp.priv, d) == nullptr) { EC_KEY_free(eckey); ossl_throw("BN_copy failed"); }
    kp.pub = EC_POINT_new(ctx.group);
    if (!kp.pub) { EC_KEY_free(eckey); ossl_throw("EC_POINT_new failed"); }
    if (EC_POINT_copy(kp.pub, Q) != 1) { EC_KEY_free(eckey); ossl_throw("EC_POINT_copy failed"); }

    EC_KEY_free(eckey);
    return kp;
}

static std::string point_hex(const OsslCtx& ctx, const EC_POINT* P) {
    char* hex = EC_POINT_point2hex(ctx.group, P, POINT_CONVERSION_UNCOMPRESSED, nullptr);
    if (!hex) ossl_throw("EC_POINT_point2hex failed");
    std::string s(hex);
    OPENSSL_free(hex);
    return s;
}

static std::string bn_hex(const BIGNUM* x) {
    char* hex = BN_bn2hex(x);
    if (!hex) ossl_throw("BN_bn2hex failed");
    std::string s(hex);
    OPENSSL_free(hex);
    return s;
}


static void mul_privates_mod_order(const OsslCtx& ctx, const std::vector<BIGNUM*>& privs, BIGNUM* out_k) {
    BN_CTX* bn = ctx.bn;
    BN_CTX_start(bn);
    BIGNUM *order = BN_CTX_get(bn);
    BIGNUM *tmp   = BN_CTX_get(bn);
    if (!tmp) { BN_CTX_end(bn); ossl_throw("BN_CTX_get failed"); }
    if (EC_GROUP_get_order(ctx.group, order, bn) != 1) { BN_CTX_end(bn); ossl_throw("EC_GROUP_get_order failed"); }

    if (BN_one(out_k) != 1) { BN_CTX_end(bn); ossl_throw("BN_one failed"); }
    for (auto* k : privs) {
        if (BN_mod_mul(tmp, out_k, k, order, bn) != 1) { BN_CTX_end(bn); ossl_throw("BN_mod_mul failed"); }
        if (BN_copy(out_k, tmp) == nullptr) { BN_CTX_end(bn); ossl_throw("BN_copy failed"); }
    }
    BN_CTX_end(bn);
}


static EC_POINT* mul_base(const OsslCtx& ctx, const BIGNUM* k) {
    BN_CTX* bn = ctx.bn;
    EC_POINT* R = EC_POINT_new(ctx.group);
    if (!R) ossl_throw("EC_POINT_new failed");
    if (EC_POINT_mul(ctx.group, R, k, nullptr, nullptr, bn) != 1) {
        EC_POINT_free(R); ossl_throw("EC_POINT_mul (base) failed");
    }
    return R;
}


static EC_POINT* mul_point(const OsslCtx& ctx, const EC_POINT* P, const BIGNUM* k) {
    BN_CTX* bn = ctx.bn;
    EC_POINT* R = EC_POINT_new(ctx.group);
    if (!R) ossl_throw("EC_POINT_new failed");
    if (EC_POINT_mul(ctx.group, R, nullptr, P, k, bn) != 1) {
        EC_POINT_free(R); ossl_throw("EC_POINT_mul (point) failed");
    }
    return R;
}


static std::vector<unsigned char> ecdh_derive_two_party_bytes() {
    std::vector<unsigned char> secret;

    auto make = []{
        EVP_PKEY_CTX* pctx = EVP_PKEY_CTX_new_id(EVP_PKEY_EC, nullptr);
        if (!pctx) ossl_throw("EVP_PKEY_CTX_new_id failed");
        if (EVP_PKEY_keygen_init(pctx) != 1) { EVP_PKEY_CTX_free(pctx); ossl_throw("keygen_init failed"); }
        if (EVP_PKEY_CTX_set_ec_paramgen_curve_nid(pctx, NID_X9_62_prime256v1) != 1) { EVP_PKEY_CTX_free(pctx); ossl_throw("set curve failed"); }
        EVP_PKEY* pkey = nullptr;
        if (EVP_PKEY_keygen(pctx, &pkey) != 1) { EVP_PKEY_CTX_free(pctx); ossl_throw("keygen failed"); }
        EVP_PKEY_CTX_free(pctx);
        return pkey;
    };

    EVP_PKEY* a = make();
    EVP_PKEY* b = make();

    EVP_PKEY_CTX* d1 = EVP_PKEY_CTX_new(a, nullptr);
    if (!d1) ossl_throw("EVP_PKEY_CTX_new failed");
    if (EVP_PKEY_derive_init(d1) != 1) { EVP_PKEY_CTX_free(d1); ossl_throw("derive_init failed"); }
    if (EVP_PKEY_derive_set_peer(d1, b) != 1) { EVP_PKEY_CTX_free(d1); ossl_throw("set_peer failed"); }
    size_t outlen=0;
    if (EVP_PKEY_derive(d1, nullptr, &outlen) != 1) { EVP_PKEY_CTX_free(d1); ossl_throw("derive(size) failed"); }
    secret.resize(outlen);
    if (EVP_PKEY_derive(d1, secret.data(), &outlen) != 1) { EVP_PKEY_CTX_free(d1); ossl_throw("derive failed"); }
    secret.resize(outlen);
    EVP_PKEY_CTX_free(d1);

    EVP_PKEY_free(a);
    EVP_PKEY_free(b);
    return secret;
}

int main() {
    try {
        OsslCtx ctx; 
        int n;
        std::cout << "Количество участников: ";
        if (!(std::cin >> n) || n < 2) {
            std::cerr << "n должно быть >= 2\n";
            return 1;
        }

        std::vector<std::unique_ptr<KeyPair>> parties;
        parties.reserve(n);
        for (int i = 0; i < n; ++i) {
            auto kp = std::make_unique<KeyPair>();
            *kp = generate_keypair(ctx);
            parties.push_back(std::move(kp));
        }

      
        std::cout << "\nПубличные ключи (uncompressed hex):\n";
        for (int i = 0; i < n; ++i) {
            std::cout << "U" << (i+1) << ": " << point_hex(ctx, parties[i]->pub) << "\n";
        }

        
        BN_CTX* bn = ctx.bn;
        BN_CTX_start(bn);
        BIGNUM* k_total = BN_CTX_get(bn);
        std::vector<BIGNUM*> privs; privs.reserve(n);
        for (int i = 0; i < n; ++i) privs.push_back(parties[i]->priv);
        mul_privates_mod_order(ctx, privs, k_total);

        EC_POINT* shared_direct = mul_base(ctx, k_total);
        std::cout << "\nОбщий секрет (прямо, prod(k_i)·G): "
                  << point_hex(ctx, shared_direct) << "\n";

       
        const EC_POINT* G = EC_GROUP_get0_generator(ctx.group);
        if (!G) { EC_POINT_free(shared_direct); BN_CTX_end(bn); ossl_throw("get generator failed"); }
        EC_POINT* step = EC_POINT_new(ctx.group);
        if (!step) { EC_POINT_free(shared_direct); BN_CTX_end(bn); ossl_throw("EC_POINT_new failed"); }
        if (EC_POINT_copy(step, G) != 1) { EC_POINT_free(step); EC_POINT_free(shared_direct); BN_CTX_end(bn); ossl_throw("EC_POINT_copy failed"); }

        for (int i = 0; i < n; ++i) {
            EC_POINT* next = mul_point(ctx, step, parties[i]->priv);
            EC_POINT_free(step);
            step = next;
        }
        std::cout << "Общий секрет (итеративно, ((((G·k1)·k2)...))): "
                  << point_hex(ctx, step) << "\n";

     
        if (EC_POINT_cmp(ctx.group, shared_direct, step, bn) != 0) {
            std::cerr << "\n[!] Несовпадение: что-то ты опять накосячил.\n";
        } else {
            std::cout << "[ok] Совпадает\n";
        }

        EC_POINT_free(shared_direct);
        EC_POINT_free(step);
        BN_CTX_end(bn);

        
        auto secret = ecdh_derive_two_party_bytes();
        std::cout << "\nДвухсторонний ECDH (EVP) секрет (" << secret.size() << " байт): ";
        for (auto b : secret) std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)b;
        std::cout << std::dec << "\n";

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Ошибка: " << ex.what() << "\n";
        return 2;
    }
}