package org.asterope.healpix;

import junit.framework.TestCase;

/**
 * @author N Kuropatkin
 *
 */
public class BitManipulationTest extends TestCase {
	/**
	 * tests bit manipulation methods
	 */
	public void testBitManipulation() {
        long mag1 = BitManipulation.magic2;
		long a = 3;
		long b = BitManipulation.swapLSBMSB(a);
		assertEquals("swapLSBMSB=" + b, 1, a/b, 1e-10);
		a = 8;
		b = BitManipulation.swapLSBMSB(a);
		assertEquals("swapLSBMSB=" + b, 2, a/b, 1e-10);
		a = 3;
		b = BitManipulation.invswapLSBMSB(a);
		assertEquals("invswapLSBMSB=" + b, -4, b, 1e-10);
		a = 8;
		b = BitManipulation.invswapLSBMSB(a);
		assertEquals("invswapLSBMSB=" + b, -5, b, 1e-10);
		
		a = 3;
		b = BitManipulation.invMSB(a);
		assertEquals("invMSB=" + b, mag1-1, b, 1e-10);
		a = 8;
		b = BitManipulation.invMSB(a);
		assertEquals("invMSB=" + b, mag1-8, b, 1e-10);
	}
	/**
	 * test Modulo
	 */
	public void testMODULO() {
		
		double a = 5.;
		double b = 3.;
		double mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = -5.0;
		b = 3.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = 5.0;
		b = -3.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = -5.0;
		b = -3.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = 8.0;
		b = 5.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = -8.0;
		b = 5.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
		a = 1.0;
		b = 4.0;
		mod = BitManipulation.MODULO(a,b);
		System.out.println("a="+a+" b="+b+" mod="+mod);
	}
}
