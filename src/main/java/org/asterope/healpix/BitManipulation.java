//
// Licenced under GPLv2, see licence.txt
// (c) K.M. Gorski, Nickolai Kuropatkin, Jan Kotek,
//
package org.asterope.healpix;

/**
 *
 *
 * bit manipulation class derived from Healpix
 * fortran90 program.
 * 
 *  @author N Kuropatkin
 */
public final class BitManipulation {
	/**
	 * magic1 odd bits set constant
	 */
	private static final long magic1 = 89478485; //  101010101010101010101010101 LSB
	/**
	 * magic2 - even bits set constant
	 */
	static final long magic2 = 178956970; //1010101010101010101010101010 MSB


	private BitManipulation() {}
	/**
	 * swaps low and high bits in the word i
	 * @param i  integer input word
	 * @return  int a word with odd and even bits interchanged
	 */
	 static public long swapLSBMSB(long i) {
		long lsb = (i & magic1);
		long msb = (i & magic2);
		return msb/2 + lsb*2;
	}
	 /**
	  * returns NOT i with even and odd bit positions interchanged
	  * @param i  int input word
	  * @return  int NOT (i with LSBMSB)
	  */
	 static public long invswapLSBMSB(long i) {

	 	long lsb = (i & magic1);
	 	long msb = (i & magic2);
	 	return ~(msb/2+lsb*2);
	 }

	 /**
	  * returns i with odd bits inverted
	  * @param i int input word
	  * @return int word with modified bits
	  */
	 static public long invMSB(long i) {
	 	return (i ^ magic2);
	 }


	 /**
	  * simulates behaviour of fortran90 MODULO function
	  * @param a  double
	  * @param b  double
	  * @return  double MODULO
	  */
	 public static double MODULO(double a, double b) {
	 	double res = 0.;
	 	long k;
	 	if (a>0.) {
	 		if (b>0.) {
	 			k = (long) (a/b);
	 			res = a - k*b;
	 			return res;
	 		}
	 		if (b<0.) {
	 			k = (long)Math.rint(a/b);
	 			res = a - k*b;
	 			return res;
	 		}
	 	}
	 	if (a<=0.) {
	 		if (b<=0.) {
	 			k = (long)(a/b);
	 			res = a - k*b;
	 			return res;
	 		}
	 		if (b>0.) {
	 			k = (long)Math.rint(a/b);
	 			res = a - k*b;
	 			return res;
	 		}
	 	}
	 	return res;
	 }
}
