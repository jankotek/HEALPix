/*
 * Created on Mar 10, 2005
 * Modified on December 2007
 *
 */
package org.asterope.healpix;

import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * test suit for PixTools
 * 
 * @author N. Kuropatkin
 * 
 */

public class PixToolsTest extends TestCase {

	/**
	 *  test MODULO function
	 */
	public void testMODULO() {
		double A = 8.;
		double B = 5.;
		double res =BitManipulation.MODULO(A, B);
		assertEquals("modulo = " + res, 3., res, 1e-10);
		A = -8.;
		B = 5.;
		res =BitManipulation.MODULO(A, B);
		assertEquals("modulo = " + res, 2., res, 1e-10);
		A = 8.;
		B = -5.;
		res =BitManipulation.MODULO(A, B);
		assertEquals("modulo = " + res, -2., res, 1e-10);
		A = -8.;
		B = -5.;
		res =BitManipulation.MODULO(A, B);
		assertEquals("modulo = " + res, -3., res, 1e-10);
		System.out.println(" test MODULO is done");
	}



	/**
	 * tests calculation of npixels from nside
	 */
	public void testNside2Npix() {
		int nside = 1;
		int npix = 0;
		npix = (int) PixTools.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 12, npix, 1e-10);
		nside = 2;
		npix = (int) PixTools.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 48, npix, 1e-10);
	}

	/**
	 * tests calculation of nsides from npixels
	 */
	public void testNpix2Nside() {
		int npix = 12;

		long nside = PixTools.Npix2Nside(npix);

		double pixSize1 = PixTools.PixRes(65536);
		long nside1 = PixTools.GetNSide(pixSize1);
		assertEquals("Nside=" + nside1, 65536, nside1, 1e-10);
		
		assertEquals("Nside=" + nside, 1, nside, 1e-10);

	}

	/**
	 * test of directional angles calculation
	 */
	public void testVec2Ang() {
		double PI = Math.PI;
		Vector3d v = new Vector3d(0.0, 1.0, 0.0);
		double[] ang_tup = { 0., 0. };
		ang_tup = PixToolsUtils.Vect2Ang(v);
		System.out.println(" Theta=" + ang_tup[0] / PI + " Phi=" + ang_tup[1]
				/ PI);
		assertEquals("Theta=" + ang_tup[0], 0.5, ang_tup[0] / PI, 1e-10);
		assertEquals("Phi=" + ang_tup[1], 0.5, ang_tup[1] / PI, 1e-10);
		v = new Vector3d(1.0, 0.0, 0.0);
		ang_tup = PixToolsUtils.Vect2Ang(v);
		assertEquals("phi=" + ang_tup[1], 0., ang_tup[1] / PI, 1e-10);
		System.out.println(" test Vect2Ang is done");
	}

	/**
	 * tests calculation of pixel from polar angles
	 * in ring schema of pixalization
	 * @throws Exception
	 */
	public void testAng2Pix() throws Exception {
		System.out.println(" Test ang2pix ___________________");
		double PI = Math.PI;
		long pix = -1;
		double theta = PI / 2. - 0.2;
		double phi = PI / 2. ; 
		long nside = 4;
		PixTools pt = new PixTools(nside);
		try {
			pix =  pt.ang2pix(theta, phi);
		} catch (Exception e) {
			e.printStackTrace();
		}
		Vector3d v = PixTools.Ang2Vec(theta,phi);
		long pix1 =  pt.vect2pix( v);
		assertEquals("pix=" + pix, pix1, pix, 1e-10);
		assertEquals("pix=" + pix, 76, pix, 1e-10);


		double[] radec = pt.pix2ang(76);
		assertEquals("theta=" + theta, theta, radec[0], 4e-2);
		assertEquals("phi=" + phi, radec[1], phi, 1e-2);
		System.out.println(" test Ang2Pix is done");
	}

	/**
	 * tests calculation of unit vector from polar angles
	 * @throws Exception
	 */
	public void testAng2Vect() throws Exception {
		System.out.println(" Start test Ang2Vect----------------");
		double PI = Math.PI;
		double theta = PI / 2.;
		double phi = PI / 2;
		Vector3d v = PixTools.Ang2Vec(theta, phi);
		System.out.println(v);
		assertEquals("x=" + v.x, 0., v.x, 1e-10);
		assertEquals("y=" + v.y, 1., v.y, 1e-10);
		assertEquals("z=" + v.z, 0., v.z, 1e-10);
		System.out.println(" test Ang2Vect is done");
	}

	/**
	 * tests calculation of ring number from z coordinate
	 * @throws Exception
	 */
	public void testRingNum() throws Exception {
		double z = 0.25;

		System.out.println("Start test RingNum !!!!!!!!!!!!!!!!!!!!");
		PixTools pt = new PixTools(1);
		int nring = (int) pt.RingNum(z);
		System.out.println("z=" + z + " ring number =" + nring);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = -0.25;
		nring = (int) pt.RingNum(z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = 0.8;
		nring = (int) pt.RingNum(z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		z = -0.8;
		nring = (int) pt.RingNum(z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		System.out.println(" test RingNum is done");
		pt = new PixTools(4);
		int pix = 3;
		Vector3d v = pt.pix2vect(pix);
		z = v.z;
		nring = (int) pt.RingNum(z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		pix = 11;
		v = pt.pix2vect(pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		pix = 23;
		v = pt.pix2vect(pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		pix = 39;
		v = pt.pix2vect(pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 4, nring, 1e-10);
		pix = 55;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 5, nring, 1e-10);
		pix = 71;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 6, nring, 1e-10);
		pix = 87;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 7, nring, 1e-10);
		pix = 103;
		v = pt.pix2vect(pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 8, nring, 1e-10);
		pix = 119;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 9, nring, 1e-10);
		pix = 135;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 10, nring, 1e-10);
		pix = 151;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 11, nring, 1e-10);
		pix = 167;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 12, nring, 1e-10);
		pix = 169;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 13, nring, 1e-10);
		pix = 180;
		v = pt.pix2vect( pix);
		z = v.z;
		nring = (int) pt.RingNum( z);
		assertEquals("z=" + z, 14, nring, 1e-10);
		System.out.println("End test RingNum");
	}




	/**
	 * tests InRing method
	 * @throws Exception
	 */
	public void testInRing() throws Exception {
		System.out.println(" Start test InRing !!!!!!!!!!!!!!!!!!!!!!!!!");
                int nside = 2;
		PixTools pt = new PixTools(nside);

		long[] ringHi = {17, 18, 19, 12, 13 };
                Arrays.sort(ringHi);

		long[] ringLow = {19, 12, 13, 14, 15 };
                    Arrays.sort(ringLow);

		double PI = Math.PI;

		int iz = 3;
		double phi = PI;
		double dphi = PI;
		LongList ring = new LongList(pt.InRing( iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " + ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		Vector3d v = new Vector3d(1., 0., 0.);
		double[] ang_tup = { 0., 0. };
		ang_tup = PixToolsUtils.Vect2Ang(v);
		phi = ang_tup[1]/PI;
		ring = new LongList(pt.InRing( iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		Vector3d v1 = new Vector3d(-1., 0., 0.);
		
		ang_tup = PixToolsUtils.Vect2Ang(v1);
		phi = ang_tup[1]/PI;
		ring = new LongList(pt.InRing( iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		phi = 1.75*PI;
		dphi = 0.5*PI;
		ring = new LongList(pt.InRing(iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					ringHi[i],  ring.get(i), 1e-10);

		}

		phi = 0.25*PI;
		dphi = 0.5*PI;

		ring =new LongList( pt.InRing( iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					ringLow[i],  ring.get(i), 1e-10);

		}			

		nside = 4;
                pt = new PixTools(nside);
		phi = 2.1598449493429825;
		iz = 8;
		dphi = 0.5890486225480867;
		//		System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = new LongList(pt.InRing( iz, phi, dphi));
		//		for (int i = 0; i<ring.size(); i++) {
		//			System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		//		}

		nside = 4;
		dphi = 0. * PI;
		iz = 8;
		phi = 2.1598449493429825;
		//		System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = new LongList(pt.InRing( iz, phi, dphi));
		//		for (int i = 0; i<ring.size(); i++) {
		//			System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		//		}
		System.out.println(" test InRing is done");
	}




	/**
	 * tests conversion from pixel number to vector
	 * @throws Exception
	 */
	public void testPix2Vect_ring() throws Exception {
		System.out.println("Start test Pix2Vect_ring !!!!!!!!!!!!!!!!!!!");
		double TWOPI = 2.0 * Math.PI;
		int nside = 2;
		int ipix = 0;
		PixTools pt = new PixTools(nside);
		Vector3d v1 = new Vector3d(0., 0., 0.);
		v1 = pt.pix2vect( ipix);
		assertEquals("v1.z = " + v1.z, 1.0, v1.z, 1e-1);

		ipix = 20;
		Vector3d v2 = new Vector3d(0., 0., 0.);
		v2 = pt.pix2vect( ipix);
		assertEquals("v2.x = " + v2.x, 1.0, v2.x, 1e-1);
		assertEquals("v2.z = " + v2.z, 0.0, v2.z, 1e-1);
		ipix = 22;
		Vector3d v3 = pt.pix2vect( ipix);
		assertEquals("v3.y = " + v3.y, 1.0, v3.y, 1e-1);
		assertEquals("v3.z = " + v3.z, 0.0, v3.z, 1e-1);
		//		System.out.println("Vector3 x="+v3.x+" y="+v3.y+" z="+v3.z);
		ipix = 95;
		nside = 4;
                pt = new PixTools(nside);
		v1 = pt.pix2vect( ipix);
		v1 = v1.normalized();
		double phi1 = Math.atan2(v1.y, v1.x);
		double[] tetphi = new double[2];
		tetphi = pt.pix2ang( ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		ipix = 26;
		nside = 4;
                pt = new PixTools(nside);
		v1 = pt.pix2vect( ipix);
		v1 = v1.normalized();
		phi1 = Math.atan2(v1.y, v1.x);
		if (phi1 < 0.)
			phi1 += TWOPI;
		tetphi = new double[2];
		tetphi = pt.pix2ang( ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		System.out.println("------------------------------------------");
		System.out.println(" test pix2vect is done");
	}




	/**
	 * tests Query_Strip method
	 * @throws Exception
	 */
	public void testQuery_Strip() throws Exception {
		System.out.println(" Start test query Strip !!!!!!!!!!!!!!!!");

		int[] pixel1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
				16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
				32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
				48, 49, 50, 51, 52, 53, 54, 55 };
		int nside = 4;
		double theta1 = 0.0;
		double theta2 = Math.PI / 4.0 + 0.2;
                PixTools pt = new PixTools(nside);
		long[] pixlist =  pt.query_strip(theta1, theta2).toArray();
		int nlist = pixlist.length;
		for (int i = 0; i < nlist; i++) {
			long ipix = pixlist[i];
			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}

		System.out.println(" test query_strip is done");

	}

	/**
	 * tests Query_Disc method
	 * @throws Exception
	 */
	public void testQuery_Disc() throws Exception {
		System.out.println(" Start test query_disc !!!!!!!!!!!!!!!!!!!!!");
                long nside = 4;
		PixTools pt = new PixTools(nside);

		long ipix = 0;
		long[] pixel1 = { 45, 46, 60, 61, 62, 77, 78,  92, 93, 94,
				 109, 110,  124, 125, 126, 141, 142 };
        Arrays.sort(pixel1);
				
		long[] pixel2 = { 24, 19, 93, 18, 17,  87, 16,86, 85,
				106,  84, 159,  81, 158, 157, 155, 156 };
        Arrays.sort(pixel2);
		long[] pixel3 = { 52, 79, 49, 78, 77,  75, 76,  74, 73, 70,
				  72, 67,  189, 66, 65, 183, 64 };
        Arrays.sort(pixel3);
		
		boolean inclusive = true;
		double radius = Math.PI / 8.0;
		Vector3d v = pt.pix2vect( 93);
		LongList pixlist;
		pixlist = new LongList(pt.query_disc(v, radius, inclusive));
		System.out.println(pixlist);
		int nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
//			System.out.println("pixel="+ipix);
			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}
		System.out.println(" End query_disk test ______________________________");
	}

	/**
	 * tests Query_Triangle method
	 * @throws Exception
	 */
	public void testQuery_Triangle() throws Exception {
                long nside = 4;
		PixTools pt = new PixTools(nside);
		long ipix = 0;
		int[] pixel1 = { 57, 58, 59, 60, 61, 62, 74, 75, 76, 77, 78, 90, 91,
				92, 93, 107, 108, 109, 123, 124, 140 };
		int[] pixel2 = { 88, 89, 90, 91,  105, 106, 107, 108, 121, 122, 123,
				138, 139, 154 };
		int[] pixel3 = { 49, 64,  80, 81, 95, 96,  112, 113, 
				127, 128, 142, 143, 144, 145 };
		int[] pixel4 = { 36, 52, 53, 67, 68, 69, 83, 84, 85, 86, 98, 99, 100,
				101, 102, 114, 115, 116, 117, 118, 119, 129, 130, 131, 132,
				133, 134, 135 };
		int[] pixel6 = { 110, 123, 124, 125, 140, 141, 156 };
		int[] pixel7 = { 53, 68, 69 };
		long pix1 = 62;
		long pix2 = 57;
		long pix3 = 140;
		System.out.println("Start test Query Triangle !!!!!!!!!!!!!!!!!!!!");
		Vector3d v11 = pt.pix2vect(pix1);

		Vector3d v22 = pt.pix2vect( pix2);

		Vector3d v33 = pt.pix2vect( pix3);

		//		System.out.println("nside="+nside+" triangle pixels "+pix1+" "+pix2+"
		// "+pix3);

		LongList pixlist = new LongList(pt.query_triangle( v11, v22, v33, false));

		int nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);

			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}
		pix1 = 92;
		pix2 = 88;
		pix3 = 154;
		v11 = pt.pix2vect( pix1);
		v22 = pt.pix2vect( pix2);
		v33 = pt.pix2vect( pix3);


		LongList pixlist1;
		pixlist1 = new LongList(pt.query_triangle( v11, v22, v33, false));

		nlist = pixlist1.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist1.get(i);
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);
		}
		pix1 = 49;
		pix2 = 142;
		pix3 = 145;
		v11 = pt.pix2vect( pix1);
		v22 = pt.pix2vect( pix2);
		v33 = pt.pix2vect( pix3);


		LongList pixlist2;
		pixlist2 = new LongList(pt.query_triangle( v11, v22, v33,  false));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix, pixel3[i], ipix, 1e-10);
		}
		pix1 = 36;
		pix2 = 129;
		pix3 = 135;
		v11 = pt.pix2vect( pix1);
		v22 = pt.pix2vect( pix2);
		v33 = pt.pix2vect( pix3);



		pixlist2 = new LongList(pt.query_triangle( v11, v22, v33,  false));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel4[i], ipix, 1e-10);
		}
		pix1 = 123;
		pix2 = 156;
		pix3 = 110;
		v11 = pt.pix2vect( pix1);
		v22 = pt.pix2vect( pix2);
		v33 = pt.pix2vect( pix3);


		pixlist2 = new LongList(pt.query_triangle( v11, v22, v33,  false));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel6[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel#="+ipix);
		}
		pix1 = 69;
		pix2 = 53;
		pix3 = 68;

		v11 = pt.pix2vect( pix1);
		v22 = pt.pix2vect( pix2);
		v33 = pt.pix2vect( pix3);


		pixlist2 = new LongList(pt.query_triangle( v11, v22, v33, false));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel7[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel#="+ipix);
		}
		System.out.println(" test query_triangle is done");

	}

	/**
	 * tests Query_Poligon method
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	public void testQuery_Polygon() throws Exception {
                long nside = 4;
		PixTools pt = new PixTools(nside);
		long ipix = 0;

		int[] result = { 51, 52, 53, 66, 67, 68, 69, 82, 83, 84, 85, 86, 98,
				99, 100, 101, 115, 116, 117 };
		int[] result1 = { 55, 70, 71, 87 };
		int[] result2 = { 137, 152, 153, 168 };
		int[] result3 = { 27, 43, 44, 58, 59, 60, 74, 75, 76, 77, 89, 90, 91,
				92, 93, 105, 106, 107, 108, 109, 110, 121, 122, 123, 124, 125,
				138, 139, 140, 141, 154, 156 };


		System.out.println("Start test query_polygon !!!!!!!!!!!!!!!!!!!!!!");
		ArrayList vlist = new ArrayList();
		Vector3d v = pt.pix2vect( 53);
		vlist.add( v);
		v = pt.pix2vect( 51);
		vlist.add( v);
		v = pt.pix2vect( 82);
		vlist.add( v);
		v = pt.pix2vect( 115);
		vlist.add( v);
		v = pt.pix2vect( 117);
		vlist.add( v);
		v = pt.pix2vect( 86);
		vlist.add( v);

		LongList pixlist = new LongList(pt.query_polygon( vlist, false));
		//		System.out.println(" List size="+pixlist.size());
		int nlist = pixlist.size();
		//		System.out.println(" Pixel list:");
		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			assertEquals("pixel = " + ipix, result[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel # "+ipix);
		}

		/* Yet another test */

		ArrayList vlist1 = new ArrayList();
		v = pt.pix2vect( 71);
		vlist1.add( v);
		v = pt.pix2vect( 55);
		vlist1.add( v);
		v = pt.pix2vect( 70);
		vlist1.add( v);
		v = pt.pix2vect( 87);
		vlist1.add( v);
		pixlist = new LongList(pt.query_polygon( vlist1, false));

		nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			//			System.out.println("i="+i+" pixel # "+ipix);
			assertEquals("pixel = " + ipix, result1[i], ipix, 1e-10);
		}

		/* Yet another test */
		ArrayList vlist2 = new ArrayList();
		v = pt.pix2vect( 153);
		vlist2.add( v);
		v = pt.pix2vect( 137);
		vlist2.add( v);
		v = pt.pix2vect( 152);
		vlist2.add( v);
		v = pt.pix2vect( 168);
		vlist2.add( v);
		pixlist = new LongList(pt.query_polygon( vlist2,  false));

		nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			assertEquals("pixel = " + ipix, result2[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel # "+ipix);
		}
		/* Yet another test */

		ArrayList vlist3 = new ArrayList();
		v = pt.pix2vect( 110);
		vlist3.add( v);
		v = pt.pix2vect( 27);
		vlist3.add( v);
		v = pt.pix2vect( 105);
		vlist3.add( v);
		v = pt.pix2vect( 154);
		vlist3.add( v);
		v = pt.pix2vect( 123);
		vlist3.add( v);
		v = pt.pix2vect( 156);
		vlist3.add( v);
		pixlist = new LongList(pt.query_polygon( vlist3,  false));

		nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			assertEquals("pixel = " + ipix, result3[i], ipix, 1e-10);
			//						System.out.println("i="+i+" pixel # "+ipix);
		}
		System.out.println(" test query_polygon is done");

	}
    /**
     * tests MaxResolution method
     */
    public void testMaxResolution() {
        System.out.println(" Start test MaxRes !!!!!!!!!!!!!!!!!!!!!");

        long nside = 1048576;
        double res = PixTools.PixRes(nside);
        System.out.println("Minimum size of the pixel side is "+res+" arcsec.");
        assertEquals("res = " + res, 0.2, res, 1e-1);
        long nsideR = PixTools.GetNSide(res);
        assertEquals("nside = " + nside, nside, nsideR, 1e-1);
        System.out.println(" End of MaxRes test _______________________");
    }
    /**
     * tests QueryDiscResolution method
     * @throws Exception
     */
    public void testQueryDiscRes() throws Exception {
        System.out.println(" Start test DiscRes !!!!!!!!!!!!!!!!!!!!!");


        long ipix = 0;
 
        boolean inclusive = false;
        double theta= Math.PI;
        double phi = Math.PI;
        double radius = Math.toRadians(0.2/3600.); //  One arcse
        long nside = PixTools.GetNSide(radius);
        PixTools pt = new PixTools(nside);
        System.out.println(" calculated nside="+nside);
        long cpix = pt.ang2pix( theta, phi);
        Vector3d vc = pt.pix2vect( cpix);
        LongList pixlist = new LongList(pt.query_disc( vc, radius, inclusive));

        int nlist = pixlist.size();
        for (int i = 0; i < nlist; i++) {
            ipix =  pixlist.get(i);
            Vector3d v = pt.pix2vect( ipix);
            double dist = v.angle(vc);
            assertTrue(dist<=2.*radius);
        }

        System.out.println(" test query disk  is done -------------------"); 
    }

    /**
     * tests GetNside method
     */
    public void testGetNside() {
        System.out.println(" Start test GetNside !!!!!!!!!!!!!!!!!!!!!");

        double pixsize = 0.3;

        long nside = PixTools.GetNSide(pixsize);
        System.out.println("Requared nside is "+nside);
        assertEquals("nside = " + nside, 1048576, nside, 1e-1);
        System.out.println(" End of GetNSide test _______________________");
    }
    

    public void testQueryCircle(){
        long nside = 512;
    	PixTools pt = new PixTools(nside);
    	


     	double angle = Math.toRadians(0.011451621372724687);
   	   	Vector3d v = new Vector3d(0.8956388362603873, -1.838600645782914E-4, 0.44478201534866);
    	
   	   	//convert vector to IPIX
    	long ipix = pt.vect2pix( v);
    	//and query circle
    	LongRangeSet r = pt.query_disc(v,angle,true);
    	
    	//now make test that IPIX is in Circle, this will fail
    	System.out.println("ipix = "+ipix);
    	System.out.println("disc: "+r);
    	assertTrue("pixel not found in disc",r.contains(ipix));

    }
}

