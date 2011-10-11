package org.asterope.healpix;

import junit.framework.TestCase;

public class PixToolsUtilsTest extends TestCase {

	/**
	 * @throws Exception
	 */
	public void testSurfaceTriangle() throws Exception {
		PixToolsVector3d v1 = new PixToolsVector3d(1.0, 0.0, 0.0);
		PixToolsVector3d v2 = new PixToolsVector3d(0.0, 1.0, 0.0);
		PixToolsVector3d v3 = new PixToolsVector3d(0.0, 0.0, 1.0);
		double res = PixToolsUtils.SurfaceTriangle(v1, v2, v3);
		System.out.println("Triangle surface is=" + res / Math.PI
				+ " steredians");
		assertEquals("Triangle surface=" + res, 0.5, res / Math.PI, 1e-10);
		System.out.println(" test of SurfaceTriangle is done");
	}
	
	 
    /**
     *  test conversion of Ra Dec to polar coordinates
     */
    public void testRaDecToPolar() {
    	System.out.println(" Start test RaDecToPolar !!!!!!!!!!!!!!!!!!!!!");
    	double [] radec = new double[2];
    	radec[0] = 312.115456;
    	radec[1] = -1.153759;
    	double[] polar = PixToolsUtils.RaDecToPolar(radec);
    	assertEquals("theta = " + polar[0], 1.5909332201194137, polar[0], 1e-10);
    	assertEquals("phi = " + polar[1], 5.447442353563491, polar[1], 1e-10);
    	System.out.println("End test RaDecToPolar__________________________");
    	
    }
    
    public void testFloor(){
    	long i1 = (long) 9.999999999d;
    	long i2 = (long)Math.floor(9.999999999d);
    	assertEquals(i1,i2);
    }
    
	/**
	 * tests intrs_intrv method
	 * @throws Exception
	 */
	public void testIntrs_Intrv() throws Exception {
		System.out.println(" test intrs_intrv !!!!!!!!!!!!!!!!!!!!!!!!!!!");
		double[] d1 = { 1.0, 9.0 };
		double[] d2 = { 3.0, 16.0 };
		double[] di;
		//		System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = PixToolsUtils.intrs_intrv(d1, d2);
		//		System.out.println("Result "+di[0]+" - "+di[1]);
		int n12 = di.length / 2;
		assertEquals("n12 = " + n12, 1, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 3.0, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 9.0, di[1], 1e-6);
		d1 = new double[] { 0.537, 4.356 };
		d2 = new double[] { 3.356, 0.8 };
		//		System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = PixToolsUtils.intrs_intrv(d1, d2);
		n12 = di.length / 2;
		assertEquals("n12 = " + n12, 2, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 0.537, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 0.8, di[1], 1e-6);
		assertEquals("di[2] = " + di[2], 3.356, di[2], 1e-6);
		assertEquals("di[1] = " + di[3], 4.356, di[3], 1e-6);

		d1 = new double[] { 2.356194490092345, 2.356194490292345 };
		d2 = new double[] { 1.251567, 4.17 };
		//		System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = PixToolsUtils.intrs_intrv(d1, d2);
		n12 = di.length / 2;
		assertEquals("n12 = " + n12, 1, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 2.35619449009, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 2.35619449029, di[1], 1e-6);

		System.out.println(" test intrs_intrv is done");
	}

}
