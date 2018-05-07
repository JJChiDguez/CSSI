#ifndef SP
#define SP

static uint8_t expn_a = 32;
static uint8_t expn_b = 20;

static uint64_t           p[] = {0xAC0E7A06FFFFFFFF, 0x12};
static uint64_t   base_mont[] = {0xD9FEFBEAD8BA0D2B, 0x4};      
// R = 2^{64 * words} mod p
static uint64_t base_mont_2[] = {0x835010E3A34C2C1C, 0x3};      
// R^2 mod p
static uint64_t          mu[] = {0xAC0E7A0700000001, 0x96F0AD1DFAEEAC43};       
// You must update mu in Arith.S: mu * p = -1 mod 2^{64 * words}
static uint64_t       zeroM[] = {0x0000000000000000, 0x0000000000000000};
static uint64_t         uno[] = {0x0000000000000001, 0x0000000000000000};       
// "Normal" One for exponentiation
static uint64_t        unoM[] = {0xD9FEFBEAD8BA0D2B, 0x4};      
// Montgomery's one
static uint64_t        negM[] = {0xD20F7E1C2745F2D4, 0xD};      
// Montgomery's minus one

static uint64_t constant_sqrt1[] = {0xAB039E81BFFFFFFF, 0x4};   
// Used for square root computation: (p - 3) / 4
static uint64_t constant_sqrt2[] = {0x56073D037FFFFFFF, 0x9};   
// Used for square root computation: (p - 1) / 2

#define WORD_N 2
#define Log2_E 35

#endif
