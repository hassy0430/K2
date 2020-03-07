/******************************************************************************/
/* Simulate KCipher-2's FSR-A register                                        */
/* In the description, alpha look up tables have had already calucrated.      */
/* I used the primitive polynominal beta(=x^8+x^7+x^6+x+1) and calucrated     */
/* alpha_0 look up table. then I simulated FSR-A.                             */
/* Reference: S. Kiyomoto, W. Shin KDDI R&D Laboratories, Inc                 */
/* "A Description of the KCipher-2 Encryption Algorithm" ISSN: 2070-1721      */
/* Takahiro, Hashimoto                                                        */
/*                                                              March 8, 2020 */
/******************************************************************************/

#include <stdio.h>
#include<inttypes.h>

/* CONSTANT VALUE */
#define PICK1		0x00000001	//for parity calculation 
#define PICK8		0x00000100	//for detect x^8

#define BETA_GF		0x000001C3	//x^8 + x^7 + x^6 + x^1 + x^0
#if 0
#define GANMA_GF	0x0000012D
#define DELTA_GF	0x0000014D
#define ZETA_GF		0x00000165
#endif

#define ALPHA_SIZE	4
#define TABLE_SIZE	256
#define GF_INDEX	8
#define FSR_A_SIZE	5

#define DEBUG 1
#define STREAM stdout

#define STR_NEWLINE "\n"

#define LOOP_A		64

#define HORIZONTAL_LINE "**************************************************"

const int32_t beta_index[ALPHA_SIZE]  = { 71,  12,   3,  24 };
#if 0
const int32_t ganma_index[ALPHA_SIZE] = { 29,  93, 156, 230 };
const int32_t delta_index[ALPHA_SIZE] = { 248, 199, 16,  34 };
const int32_t zeta_index[ALPHA_SIZE]  = { 16,  56, 253, 157 };
#endif

/* INLINE MACRO FUNCTION */
inline uint32_t parity8(x) {
	uint32_t p;
	p = (x);
	p ^= (p >> 4);
	p ^= (p >> 2);
	p ^= (p >> 1);
	return p & PICK1;
}

inline uint32_t parity32(uint32_t x) {
	uint32_t p;
	p = (x);
	p ^= (p >> 16);
	p ^= (p >>  8);
	p ^= (p >>  4);
	p ^= (p >>  2);
	p ^= (p >>  1);
	return p & PICK1;
}

/* Multiplication in Galois Field */
inline uint32_t multi_GF(uint32_t a, uint32_t num, uint32_t gf, uint32_t gf_idx)
{
	uint32_t t = 0;
	int debug = 0;
	int pos;
	for (pos = 7; pos >= 0; pos--) {
		if ((num) & (PICK1 << pos))break;
	}
	for (int i = pos; i >= 0; i--) {
		t = (t << 1) ^ (((num) & (PICK1 << i)) ? (a) : 0);
		if (t & (PICK1 << (gf_idx)))
			t ^= (gf);
	}
	return t;
}

/* Make coefficients of a monic polynomial  */
void make_coefficients(uint32_t multi[], uint32_t index[], uint32_t gf)
{
	for (int i = 0; i < ALPHA_SIZE; i++) {
		for (int j = 0; j < (int)index[i]; j++) {
			multi[i] <<= 1;
			if (multi[i] & PICK8)
				multi[i] ^= gf | !(parity8(multi[i] ^ gf));
		}
	}
}

/* Make alpha look up table */
void make_alpha_table(uint32_t alpha[], uint32_t multi[], uint32_t gf)
{
	uint32_t temp_u32;
	for (int i = 0; i < TABLE_SIZE; i++) {
		for (int j = 0; j < ALPHA_SIZE; j++) {
			temp_u32 = multi_GF(multi[j], (uint32_t)i, gf, GF_INDEX);
			alpha[i] |= temp_u32 << (j * 8);
		}
	}
}

/* Update FSR-A registers */
void update_fsr(uint32_t fsr[], uint32_t alpha[])
{
	uint32_t fb;
	fb = (fsr[0] << 8) ^ alpha[(fsr[0] >> 24)] ^ fsr[3];
	fsr[0] = fsr[1];
	fsr[1] = fsr[2];
	fsr[2] = fsr[3];
	fsr[3] = fsr[4];
	fsr[4] = fb;
}

/* Print FSR-A registers */
void printf_fsr(uint32_t fsr[], int loop, FILE* stream)
{
	fprintf_s(stream, "%s%s", HORIZONTAL_LINE, STR_NEWLINE);
	fprintf_s(stream, "loop:%2d%s", loop, STR_NEWLINE);
	for (int i = 0; i < FSR_A_SIZE; i++) {
		fprintf_s(stream, "FSR-A[%d]:%08X%s", i, fsr[i], STR_NEWLINE);
	}
}

int main()
{
	uint32_t beta_multi[ALPHA_SIZE] = { 1, 1, 1, 1 };
	uint32_t alpha_0[TABLE_SIZE] = { 0 };
	uint32_t fsr_a[FSR_A_SIZE] = { 0xBE3CA984, 0x974E6719, 0x86916EFF, 0xF52DACF9, 0x960329B5 };
	int loop = LOOP_A;

	int debug = DEBUG;
	FILE* stream = STREAM;
	debug = 0;
	/*  Make coefficients of a monic polynomial */
	make_coefficients(beta_multi, beta_index, BETA_GF);
	for (int i = 0; i < ALPHA_SIZE; i++) {
		if (debug)fprintf_s(stream, "beta^%d = %02X%s", beta_index[i], beta_multi[i], STR_NEWLINE);
	}
	debug = 1;
	/* Make alpha look up table */
	if (debug)fprintf_s(stream, "alpha_0[256]={%s", STR_NEWLINE);
	make_alpha_table(alpha_0, beta_multi, BETA_GF);
	for (int i = 0; i < TABLE_SIZE; i++) {
		if (debug)fprintf_s(stream, "0x%08X%s%s", alpha_0[i], (i < TABLE_SIZE - 1) ? "," : "", (i + 1) % 7 ? "" : STR_NEWLINE);
	}
	if (debug)fprintf_s(stream, "};%s", STR_NEWLINE);
	
	/* FSR-A */
	if (debug)printf_fsr(fsr_a, 0, stream);
	for (int i = 0; i < loop; i++) {
		update_fsr(fsr_a, alpha_0);
		if (debug)printf_fsr(fsr_a, i + 1, stream);
	}
	if(debug)fprintf_s(stream, "%s%s", HORIZONTAL_LINE, STR_NEWLINE);

	//system("pause");

	return 0;
}
