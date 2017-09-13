/* A version of Marsaglia-MultiCarry */

static unsigned int I1 = 1234, I2 = 5678;

void set_seed(unsigned int i1, unsigned int i2) {
    I1 = i1;
    I2 = i2;
}

void get_seed(unsigned int *i1, unsigned int *i2) {
    *i1 = I1;
    *i2 = I2;
}

double unif_rand(void) {
    I1 = 36969 * (I1 & 0177777) + (I1 >> 16);
    I2 = 18000 * (I2 & 0177777) + (I2 >> 16);
    return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}
