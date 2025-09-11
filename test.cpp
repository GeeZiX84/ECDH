#include <openssl/ssl.h>
#include <openssl/err.h>
int main() {
    SSL_library_init();
    return 0;
}
