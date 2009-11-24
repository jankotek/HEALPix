/*
 * Created on Mar 10, 2005
 * Modified on December 2007
 *
 */
package org.asterope.healpix;

import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * test suit for PixToolsNested
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
		PixToolsNested pt = new PixToolsNested();
		npix = (int) pt.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 12, npix, 1e-10);
		nside = 2;
		npix = (int) pt.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 48, npix, 1e-10);
	}

	/**
	 * tests calculation of nsides from npixels
	 */
	public void testNpix2Nside() {
		int npix = 12;
		int nside = 0;
		PixToolsNested pt = new PixToolsNested();
		nside = (int) pt.Npix2Nside(npix);

		double pixSize1 = pt.PixRes(65536);
		long nside1 = pt.GetNSide(pixSize1);
		assertEquals("Nside=" + nside1, 65536, nside1, 1e-10);
		
		assertEquals("Nside=" + nside, 1, nside, 1e-10);

	}

	/**
	 * test of directional angles calculation
	 */
	public void testVec2Ang() {
		double PI = Math.PI;
		PixToolsVector3d v = new PixToolsVector3d(0.0, 1.0, 0.0);
		double[] ang_tup = { 0., 0. };
		ang_tup = PixToolsUtils.Vect2Ang(v);
		System.out.println(" Theta=" + ang_tup[0] / PI + " Phi=" + ang_tup[1]
				/ PI);
		assertEquals("Theta=" + ang_tup[0], 0.5, ang_tup[0] / PI, 1e-10);
		assertEquals("Phi=" + ang_tup[1], 0.5, ang_tup[1] / PI, 1e-10);
		v = new PixToolsVector3d(1.0, 0.0, 0.0);
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
		System.out.println(" Test ang2pix_ring ___________________");
		double PI = Math.PI;
		long pix = -1;
		double theta = PI / 2. - 0.2;
		double phi = PI / 2. ; 
		long nside = 4;
		PixToolsNested pt = new PixToolsNested();
		try {
			pix =  pt.ang2pix_ring(nside,theta, phi);
		} catch (Exception e) {
			e.printStackTrace();
		}
		PixToolsVector3d v = pt.Ang2Vec(theta,phi);
		long pix1 =  pt.vect2pix_ring(nside,v);
		assertEquals("pix=" + pix, pix1, pix, 1e-10);
		assertEquals("pix=" + pix, 76, pix, 1e-10);

		long pix2 =  pt.ang2pix_nest(nside,theta,phi);
		long pix3 =  pt.vect2pix_nest(nside,v);
		assertEquals("pix2=" + pix2, pix3, pix2, 1e-10);
		assertEquals("pix2=" + pix2, 92, pix2, 1e-10);


		double[] radec = pt.pix2ang_ring(nside,76);
		assertEquals("theta=" + theta, theta, radec[0], 4e-2);
		assertEquals("phi=" + phi, radec[1], phi, 1e-2);
		double[] radec1 = pt.pix2ang_nest(nside,92);
		System.out.println("theta="+radec1[0]+" phi="+radec1[1]);
		assertEquals("theta=" + theta, theta, radec1[0], 4e-2);
		assertEquals("phi=" + phi, radec1[1], phi, 1e-2);	
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
		PixToolsNested pt = new PixToolsNested();
		PixToolsVector3d v = pt.Ang2Vec(theta, phi);
		System.out.println("Vector x=" + v.x + " y=" + v.y + " z=" + v.z);
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
		int nside = 1;
		System.out.println("Start test RingNum !!!!!!!!!!!!!!!!!!!!");
		PixToolsNested pt = new PixToolsNested();
		int nring = (int) pt.RingNum(nside, z);
		System.out.println("z=" + z + " ring number =" + nring);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = -0.25;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = 0.8;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		z = -0.8;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		System.out.println(" test RingNum is done");
		nside = 4;
		int pix = 3;
		PixToolsVector3d v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		pix = 11;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		pix = 23;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		pix = 39;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 4, nring, 1e-10);
		pix = 55;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 5, nring, 1e-10);
		pix = 71;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 6, nring, 1e-10);
		pix = 87;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 7, nring, 1e-10);
		pix = 103;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 8, nring, 1e-10);
		pix = 119;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 9, nring, 1e-10);
		pix = 135;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 10, nring, 1e-10);
		pix = 151;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 11, nring, 1e-10);
		pix = 167;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 12, nring, 1e-10);
		pix = 169;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 13, nring, 1e-10);
		pix = 180;
		v = pt.pix2vect_ring(nside, pix);
		z = v.z;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 14, nring, 1e-10);
		System.out.println("End test RingNum");
	}

	/**
	 * tests conversion from nest schema pixel to ring schema pixel
	 * @throws Exception
	 */
	public void testNest2Ring() throws Exception {
		PixToolsNested pt = new PixToolsNested();
		int ipnest = 3;
		int nside = 2;
		int ipring = (int) pt.nest2ring(nside, ipnest);
		assertEquals("ipring=" + ipring, 0, ipring, 1e-10);
		ipnest = 0;
		nside = 2;
		ipring = (int) pt.nest2ring(nside, ipnest);
		assertEquals("ipring=" + ipring, 13, ipring, 1e-10);
		ipnest = 18;
		nside = 2;
		ipring = (int) pt.nest2ring(nside, ipnest);
		assertEquals("ipring=" + ipring, 27, ipring, 1e-10);
		ipnest = 23;
		nside = 2;
		ipring = (int) pt.nest2ring(nside, ipnest);
		assertEquals("ipring=" + ipring, 14, ipring, 1e-10);
		ipnest = 5;
		nside = 4;
		ipring = (int) pt.nest2ring(nside, ipnest);
		assertEquals("ipring = " + ipring, 27, ipring, 1e-10);
		System.out.println(" test Nest2Ring is done");
	}

	/**
	 * tests conversion from ring schema pixel to nest schema pixel
	 * @throws Exception
	 */
	public void testRing2Nest() throws Exception {
		PixToolsNested pt = new PixToolsNested();
		System.out.println(" start test Ring2Nest !!!!!!!!!!!!!!!!!!!!!!");
		int ipring = 0;
		int nside = 2;

		int ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest=" + ipnest, 3, ipnest, 1e-10);
		ipring = 13;
		nside = 2;
		ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest=" + ipnest, 0, ipnest, 1e-10);
		ipring = 27;
		nside = 2;
		ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest=" + ipnest, 18, ipnest, 1e-10);
		ipring = 14;
		nside = 2;
		ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest=" + ipnest, 23, ipnest, 1e-10);
		ipring = 27;
		nside = 4;
		ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest = " + ipnest, 5, ipnest, 1e-10);
		ipring = 83;
		nside = 4;
		ipnest = (int) pt.ring2nest(nside, ipring);
		assertEquals("ipnest = " + ipnest, 123, ipnest, 1e-10);
		System.out.println(" test Ring2Nest is done");
	}

	/**
	 * tests Next_In_Line method for the nest schema
	 * @throws Exception
	 */
	public void testNext_In_Line_Nest() throws Exception {
		PixToolsNested pt = new PixToolsNested();
		int ipix = 0;
		int nside = 2;
		int ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 23, ipnext, 1e-10);
		ipix = 1;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 6, ipnext, 1e-10);
		ipix = 4;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 27, ipnext, 1e-10);
		ipix = 27;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 8, ipnext, 1e-10);
		ipix = 12;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 19, ipnext, 1e-10);
		ipix = 118;
		nside = 4;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext = " + ipnext, 117, ipnext, 1e-10);
		System.out.println(" test next_in_line_nest is done");
	}

	/**
	 * tests InRing method
	 * @throws Exception
	 */
	public void testInRing() throws Exception {
		System.out.println(" Start test InRing !!!!!!!!!!!!!!!!!!!!!!!!!");
		PixToolsNested pt = new PixToolsNested();
		int[] nestComp = { 0,4,8,12,19,23,27,31 };
		long[] ringHi = {17, 18, 19, 12, 13 };
        Arrays.sort(ringHi);
		long[] nestHi = { 0,  8,12,  19 ,  31  };
		long[] ringLow = {19, 12, 13, 14, 15 };
        Arrays.sort(ringLow);
		long[] nestLow = {0,4,12,19,23  };
		double PI = Math.PI;
		int nside = 2;
		int iz = 3;
		double phi = PI;
		double dphi = PI;
		LongList ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " + ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		PixToolsVector3d v = new PixToolsVector3d(1., 0., 0.);
		double[] ang_tup = { 0., 0. };
		ang_tup = PixToolsUtils.Vect2Ang(v);
		phi = ang_tup[1]/PI;
		ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		PixToolsVector3d v1 = new PixToolsVector3d(-1., 0., 0.);
		
		ang_tup = PixToolsUtils.Vect2Ang(v1);
		phi = ang_tup[1]/PI;
		ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					i + 12,  ring.get(i), 1e-10);
		}
		phi = 1.75*PI;
		dphi = 0.5*PI;
		ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					ringHi[i],  ring.get(i), 1e-10);

		}

		phi = 1.75*PI;
		dphi = 0.5*PI;

		ring = new LongList(pt.InRing_nested(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					nestHi[i],  ring.get(i), 1e-10);

		}	
		phi = 0.25*PI;
		dphi = 0.5*PI;

		ring =new LongList( pt.InRing(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					ringLow[i],  ring.get(i), 1e-10);

		}			
		phi = 0.25*PI;
		dphi = 0.5*PI;
		
		ring = new LongList(pt.InRing_nested(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					nestLow[i],  ring.get(i), 1e-10);

		}	
		
		dphi = PI;
		ring = new LongList(pt.InRing_nested(nside, iz, phi, dphi));
		for (int i = 0; i < ring.size(); i++) {
			assertEquals("ipnext = " +  ring.get(i),
					nestComp[i],  ring.get(i), 1e-10);
		}

		nside = 4;
		phi = 2.1598449493429825;
		iz = 8;
		dphi = 0.5890486225480867;
		//		System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		//		for (int i = 0; i<ring.size(); i++) {
		//			System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		//		}

		nside = 4;
		dphi = 0. * PI;
		iz = 8;
		phi = 2.1598449493429825;
		//		System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = new LongList(pt.InRing(nside, iz, phi, dphi));
		//		for (int i = 0; i<ring.size(); i++) {
		//			System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		//		}
		System.out.println(" test InRing is done");
	}

	/**
	 * tests Neighbour's method for nest schema of the pixelization 
	 * @throws Exception
	 */
	public void testNeighbours_Nest() throws Exception {
		System.out.println(" Start test Neighbours_Nest !!!!!!!!!!!!!!!!!");
		PixToolsNested pt = new PixToolsNested();
		long nside = 2;
		long ipix = 25;
		long[] test17 = { 34, 16, 18, 19, 2, 0, 22, 35 };
		long[] test32 = { 40, 44, 45, 34, 35, 33, 38, 36 };
		long[] test3 = { 0, 2, 13, 15, 11, 7, 6, 1 };
		long[] test25 = { 42, 24, 26, 27, 10, 8, 30, 43 };
		//
		LongList npixList = new LongList();
		npixList = pt.neighbours_nest(nside, ipix);
		for (int i = 0; i < npixList.size(); i++) {
			assertEquals("ip = " +  npixList.get( i),
					test25[ i],  npixList.get(  i), 1e-10);
		}
		ipix = 17;

		npixList = pt.neighbours_nest(nside, ipix);
		for (int i = 0; i < npixList.size(); i++) {
			assertEquals("ip = " +  npixList.get(i),
					test17[i],  npixList.get(i), 1e-10);
		}
		ipix = 32;

		npixList = pt.neighbours_nest(nside, ipix);
		for (int i = 0; i < npixList.size(); i++) {
			assertEquals("ip = " +  npixList.get(i),
					test32[i],  npixList.get(i), 1e-10);
		}
		ipix = 3;

		npixList = pt.neighbours_nest(nside, ipix);
		for (int i = 0; i < npixList.size(); i++) {
			assertEquals("ip = " +  npixList.get(i),
					test3[i],  npixList.get(i), 1e-10);
		}
		System.out.println(" test NeighboursNest is done");
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
		PixToolsNested pt = new PixToolsNested();
		PixToolsVector3d v1 = new PixToolsVector3d(0., 0., 0.);
		v1 = pt.pix2vect_ring(nside, ipix);
		assertEquals("v1.z = " + v1.z, 1.0, v1.z, 1e-1);

		ipix = 20;
		PixToolsVector3d v2 = new PixToolsVector3d(0., 0., 0.);
		v2 = pt.pix2vect_ring(nside, ipix);
		assertEquals("v2.x = " + v2.x, 1.0, v2.x, 1e-1);
		assertEquals("v2.z = " + v2.z, 0.0, v2.z, 1e-1);
		ipix = 22;
		PixToolsVector3d v3 = pt.pix2vect_ring(nside, ipix);
		assertEquals("v3.y = " + v3.y, 1.0, v3.y, 1e-1);
		assertEquals("v3.z = " + v3.z, 0.0, v3.z, 1e-1);
		//		System.out.println("Vector3 x="+v3.x+" y="+v3.y+" z="+v3.z);
		ipix = 95;
		nside = 4;
		v1 = pt.pix2vect_ring(nside, ipix);
		v1 = v1.normalized();
		double phi1 = Math.atan2(v1.y, v1.x);
		double[] tetphi = new double[2];
		tetphi = pt.pix2ang_ring(nside, ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		ipix = 26;
		nside = 4;
		v1 = pt.pix2vect_ring(nside, ipix);
		v1 = v1.normalized();
		phi1 = Math.atan2(v1.y, v1.x);
		if (phi1 < 0.)
			phi1 += TWOPI;
		tetphi = new double[2];
		tetphi = pt.pix2ang_ring(nside, ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		System.out.println("------------------------------------------");
		System.out.println(" test pix2vect_ring is done");
	}

	/**
	 * tests conversion from pixel number to vector
	 * @throws Exception
	 */
	public void testPix2Vect_nest() throws Exception {
		double TWOPI = 2.0 * Math.PI;
		int nside = 2;
		int ipix = 3;
		System.out.println(" Start test Pix2Vect_nest !!!!!!!!!!!!!!");
		PixToolsNested pt = new PixToolsNested();
		PixToolsVector3d v1 = new PixToolsVector3d(0., 0., 0.);
		v1 = pt.pix2vect_nest(nside, ipix);
		assertEquals("v1.z = " + v1.z, 1.0, v1.z, 1e-1);
		ipix = 17;
		PixToolsVector3d v2 = new PixToolsVector3d(0., 0., 0.);
		v2 = pt.pix2vect_nest(nside, ipix);
		assertEquals("v2.x = " + v2.x, 1.0, v2.x, 1e-1);
		assertEquals("v2.z = " + v2.z, 0.0, v2.z, 1e-1);

		ipix = 21;
		PixToolsVector3d v3 = pt.pix2vect_nest(nside, ipix);
		assertEquals("v3.y = " + v3.y, 1.0, v3.y, 1e-1);
		assertEquals("v3.z = " + v3.z, 0.0, v3.z, 1e-1);
		nside = 4;
		ipix = 105;
		v1 = pt.pix2vect_nest(nside, ipix);
		v1 = v1.normalized();
		double phi1 = Math.atan2(v1.y, v1.x);
		if (phi1 < 0.)
			phi1 += TWOPI;
		double[] tetphi = new double[2];
		tetphi = pt.pix2ang_nest(nside, ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		nside = 4;
		ipix = 74;
		v1 = pt.pix2vect_nest(nside, ipix);
		v1 = v1.normalized();
		phi1 = Math.atan2(v1.y, v1.x);
		if (phi1 < 0.)
			phi1 += TWOPI;
		tetphi = new double[2];
		tetphi = pt.pix2ang_nest(nside, ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);

		System.out.println(" test pix2vect_nest is done");
		System.out.println("-------------------------------------------");
	}

	/**
	 * tests conversion from vector to pixel number
	 * @throws Exception
	 */
	public void testVect2Pix_ring() throws Exception {
		System.out.println("Start test Vect2Pix_ring !!!!!!!!!!!!!!!!!!!");

		long nside = 4;
		long ipix = 83;
		long respix = 0;
		PixToolsNested pt = new PixToolsNested();
		PixToolsVector3d v1 = new PixToolsVector3d(0., 0., 0.);
		v1 = pt.pix2vect_ring(nside, ipix);
		respix = (int) pt.vect2pix_ring(nside, v1);
		assertEquals("respix = " + respix, 83, respix, 1e-10);
		// Hi resolution test
		long nside1 = 1 << 20;
		long maxpix= pt.Nside2Npix(nside1);
		System.out.println("nside="+nside1+" maxpix="+maxpix);
		PixToolsVector3d v2 = new PixToolsVector3d( -0.704, 0.580, 0.408);
		respix = pt.vect2pix_ring(nside1, v2);  // ring pixel
		long respixN = pt.ring2nest(nside1, respix); // convert to nest
		long respixNC = pt.vect2pix_nest(nside1,v2); // nest pixel from the same vector
		long respixR = pt.nest2ring(nside1, respixN); // convert pixel 
		System.out.println(" orig="+respix+" doubleT="+respixR+" nest="+respixN+" correct nest="+respixNC);
		assertEquals("ringpix = " + respix, respix, respixR);
		assertEquals("nestpix = " + respixNC, respixNC, respixN);
		System.out.println("------------------------------------------");
		System.out.println(" test vect2pix_ring is done");
	}

	/**
	 * tests conversion from vector to pixel number
	 * @throws Exception
	 */
	public void testVect2Pix_nest() throws Exception {
		System.out.println("Start test Vect2Pix_nest !!!!!!!!!!!!!!!!!!!");
		long nside = 4;
		long ipix = 83;
		long respix = 0;
		PixToolsNested pt = new PixToolsNested();
		PixToolsVector3d v1 = new PixToolsVector3d(0., 0., 0.);
		v1 = pt.pix2vect_ring(nside, ipix);
		respix = (int) pt.vect2pix_nest(nside, v1);
		assertEquals("respix = " + respix, 123, respix, 1e-10);
		//
		long nside1 = 1 << 20;
		long maxpix= pt.Nside2Npix(nside1);
        System.out.println("nside="+nside1+" maxpix="+maxpix);
		PixToolsVector3d v2 = new PixToolsVector3d( -0.704, 0.580, 0.408);
		respix = pt.vect2pix_nest(nside1, v2);
		long respixRC = pt.vect2pix_ring(nside1, v2);
		long respixR = pt.nest2ring(nside1, respix);
		long respixN = pt.ring2nest(nside1, respixRC);
		System.out.println(" orig="+respix+" doubleT="+respixN+" ring="+respixR+" correct ring="+respixRC);
		assertEquals("ringpix = " + respixRC, respixRC, respixR);
        assertEquals("nestpix = " + respix, respix, respixN);
		System.out.println("------------------------------------------");
		System.out.println(" test vect2pix_nest is done");
	}

	/**
	 * tests Query_Strip method
	 * @throws Exception
	 */
	public void testQuery_Strip() throws Exception {
		System.out.println(" Start test query Strip !!!!!!!!!!!!!!!!");
		PixToolsNested pt = new PixToolsNested();
		int[] pixel1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
				16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
				32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
				48, 49, 50, 51, 52, 53, 54, 55 };
		int nside = 4;
		double theta1 = 0.0;
		double theta2 = Math.PI / 4.0 + 0.2;
		long[] pixlist =  pt.query_strip(nside, theta1, theta2).toArray();
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
		PixToolsNested pt = new PixToolsNested();
		long nside = 4;
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
		PixToolsVector3d v = pt.pix2vect_ring(nside, 93);
		LongList pixlist;
		pixlist = new LongList(pt.query_disc(nside, v, radius, inclusive));
		System.out.println(pixlist);
		int nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
//			System.out.println("pixel="+ipix);
			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}
		radius = Math.PI/2.0; 

		v = new PixToolsVector3d(1.0,0.0,0.0);
//		LongList pixlistT = new LongList(pt.query_disc(2, v, radius, inclusive));

		radius = Math.PI / 8.0;
		v = pt.pix2vect_ring(nside, 93);
		LongList pixlist2 = new LongList(pt.query_disc_nested(nside, v, radius,  inclusive));
		System.out.println(pixlist2);

		int nlist2 = pixlist2.size();
        //pixlist2 is sorted from LongRangeSet, sort pixel2 as well
        assertEquals("npix="+nlist, nlist, nlist2, 1e-10);
		for (int i = 0; i < nlist2; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);

		}

		v = pt.pix2vect_ring(nside, 103);
		inclusive = true;
		LongList pixlist3 = new LongList(pt.query_disc_nested(nside, v, radius, inclusive));
		nlist = pixlist3.size();
		assertEquals("npix="+nlist, nlist, pixel3.length, 1e-10);
		for (int i = 0; i < pixlist3.size(); i++) {
			ipix =  pixlist3.get(i);
//						System.out.println("ipix="+ ipix);
			assertEquals("pixel = " + ipix, pixel3[i], ipix, 1e-10);
		}

		for (int i=0; i<pixel1.length; i++) {
			long ipixR = pixel1[i];	
			long ipixT = pt.ring2nest(nside, ipixR);
            assertTrue("pixel="+ipixT, new LongList(pixel2).contains(ipixT));
			//assertEquals("pixel="+ipixT, ipixT, pixel2[i], 1e-10);
			long ipixN = pt.nest2ring(nside, ipixT);
			assertEquals("pixel="+ipixN, ipixN, pixel1[i], 1e-10);
		}
		System.out.println(" End query_disk test ______________________________");
	}

	/**
	 * tests Query_Triangle method
	 * @throws Exception
	 */
	public void testQuery_Triangle() throws Exception {
		PixToolsNested pt = new PixToolsNested();
		long nside = 4;
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
		int[] pixel5 = { 58, 127, 56, 126, 125, 50, 123, 124, 119, 48, 122,
				121, 118, 117, 74, 175, 120, 115, 116, 191, 72, 174, 173, 114,
				113, 190, 189, 66 };
		int[] pixel6 = { 110, 123, 124, 125, 140, 141, 156 };
		int[] pixel7 = { 53, 68, 69 };
		long pix1 = 62;
		long pix2 = 57;
		long pix3 = 140;
		System.out.println("Start test Query Triangle !!!!!!!!!!!!!!!!!!!!");
		PixToolsVector3d v11 = pt.pix2vect_ring(nside, pix1);

		PixToolsVector3d v22 = pt.pix2vect_ring(nside, pix2);

		PixToolsVector3d v33 = pt.pix2vect_ring(nside, pix3);

		//		System.out.println("nside="+nside+" triangle pixels "+pix1+" "+pix2+"
		// "+pix3);
		int inclusive = 0;

		LongList pixlist = new LongList(pt.query_triangle(nside, v11, v22, v33, inclusive));

		int nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);

			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}
		pix1 = 92;
		pix2 = 88;
		pix3 = 154;
		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);

		inclusive = 0;
		LongList pixlist1;
		pixlist1 = new LongList(pt.query_triangle(nside, v11, v22, v33, inclusive));

		nlist = pixlist1.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist1.get(i);
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);
		}
		pix1 = 49;
		pix2 = 142;
		pix3 = 145;
		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);

		inclusive = 0;
		LongList pixlist2;
		pixlist2 = new LongList(pt.query_triangle(nside, v11, v22, v33,  inclusive));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix, pixel3[i], ipix, 1e-10);
		}
		pix1 = 36;
		pix2 = 129;
		pix3 = 135;
		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);

		inclusive = 0;

		pixlist2 = new LongList(pt.query_triangle(nside, v11, v22, v33,  inclusive));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel4[i], ipix, 1e-10);
		}
		pix1 = 36;
		pix2 = 129;
		pix3 = 135;
		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);
		inclusive = 0;

		pixlist2 = new LongList(pt.query_triangle_nested(nside, v11, v22, v33, inclusive));

		nlist = pixlist2.size();

        //pixlist2 is sorted from LongRangeSet, sort pixel5 as well
        Arrays.sort(pixel5);
		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
//            System.out.println("ipix="+ipix);
			assertEquals("pixel = " + ipix, pixel5[i], ipix, 1e-10);
		}
		pix1 = 123;
		pix2 = 156;
		pix3 = 110;
		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);
		inclusive = 0;

		pixlist2 = new LongList(pt.query_triangle(nside, v11, v22, v33,  inclusive));

		nlist = pixlist2.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist2.get(i);
			assertEquals("pixel = " + ipix, pixel6[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel#="+ipix);
		}
		pix1 = 69;
		pix2 = 53;
		pix3 = 68;

		v11 = pt.pix2vect_ring(nside, pix1);
		v22 = pt.pix2vect_ring(nside, pix2);
		v33 = pt.pix2vect_ring(nside, pix3);
		inclusive = 0;

		pixlist2 = new LongList(pt.query_triangle(nside, v11, v22, v33, inclusive));

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
		PixToolsNested pt = new PixToolsNested();
		long nside = 4;
		long ipix = 0;
		int inclusive = 0;
		int[] result = { 51, 52, 53, 66, 67, 68, 69, 82, 83, 84, 85, 86, 98,
				99, 100, 101, 115, 116, 117 };
		int[] result1 = { 55, 70, 71, 87 };
		int[] result2 = { 137, 152, 153, 168 };
		int[] result3 = { 27, 43, 44, 58, 59, 60, 74, 75, 76, 77, 89, 90, 91,
				92, 93, 105, 106, 107, 108, 109, 110, 121, 122, 123, 124, 125,
				138, 139, 140, 141, 154, 156 };


		System.out.println("Start test query_polygon !!!!!!!!!!!!!!!!!!!!!!");
		ArrayList vlist = new ArrayList();
		PixToolsVector3d v = pt.pix2vect_ring(nside, 53);
		vlist.add( v);
		v = pt.pix2vect_ring(nside, 51);
		vlist.add( v);
		v = pt.pix2vect_ring(nside, 82);
		vlist.add( v);
		v = pt.pix2vect_ring(nside, 115);
		vlist.add( v);
		v = pt.pix2vect_ring(nside, 117);
		vlist.add( v);
		v = pt.pix2vect_ring(nside, 86);
		vlist.add( v);

		LongList pixlist = new LongList(pt.query_polygon(nside, vlist, inclusive));
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
		v = pt.pix2vect_ring(nside, 71);
		vlist1.add( v);
		v = pt.pix2vect_ring(nside, 55);
		vlist1.add( v);
		v = pt.pix2vect_ring(nside, 70);
		vlist1.add( v);
		v = pt.pix2vect_ring(nside, 87);
		vlist1.add( v);
		pixlist = new LongList(pt.query_polygon(nside, vlist1, inclusive));

		nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			//			System.out.println("i="+i+" pixel # "+ipix);
			assertEquals("pixel = " + ipix, result1[i], ipix, 1e-10);
		}

		/* Yet another test */
		ArrayList vlist2 = new ArrayList();
		v = pt.pix2vect_ring(nside, 153);
		vlist2.add( v);
		v = pt.pix2vect_ring(nside, 137);
		vlist2.add( v);
		v = pt.pix2vect_ring(nside, 152);
		vlist2.add( v);
		v = pt.pix2vect_ring(nside, 168);
		vlist2.add( v);
		pixlist = new LongList(pt.query_polygon(nside, vlist2,  inclusive));

		nlist = pixlist.size();

		for (int i = 0; i < nlist; i++) {
			ipix =  pixlist.get(i);
			assertEquals("pixel = " + ipix, result2[i], ipix, 1e-10);
			//			System.out.println("i="+i+" pixel # "+ipix);
		}
		/* Yet another test */

		ArrayList vlist3 = new ArrayList();
		v = pt.pix2vect_ring(nside, 110);
		vlist3.add( v);
		v = pt.pix2vect_ring(nside, 27);
		vlist3.add( v);
		v = pt.pix2vect_ring(nside, 105);
		vlist3.add( v);
		v = pt.pix2vect_ring(nside, 154);
		vlist3.add( v);
		v = pt.pix2vect_ring(nside, 123);
		vlist3.add( v);
		v = pt.pix2vect_ring(nside, 156);
		vlist3.add( v);
		pixlist = new LongList(pt.query_polygon(nside, vlist3,  inclusive));

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
        PixToolsNested pt = new PixToolsNested();
        double res = pt.PixRes(nside);
        System.out.println("Minimum size of the pixel side is "+res+" arcsec.");
        assertEquals("res = " + res, 0.2, res, 1e-1);
        long nsideR = pt.GetNSide(res);
        assertEquals("nside = " + nside, nside, nsideR, 1e-1);
        System.out.println(" End of MaxRes test _______________________");
    }
    /**
     * tests QueryDiscResolution method
     * @throws Exception
     */
    public void testQueryDiscRes() throws Exception {
        System.out.println(" Start test DiscRes !!!!!!!!!!!!!!!!!!!!!");
        PixToolsNested pt = new PixToolsNested();

        long ipix = 0;
 
        boolean inclusive = false;
        double theta= Math.PI;
        double phi = Math.PI;
        double radius = Math.toRadians(0.2/3600.); //  One arcse
        long nside = pt.GetNSide(radius);
        System.out.println(" calculated nside="+nside);
        long cpix = pt.ang2pix_ring(nside,theta,phi);
        PixToolsVector3d vc = pt.pix2vect_ring(nside, cpix);
        LongList pixlist = new LongList(pt.query_disc(nside, vc, radius, inclusive));

        int nlist = pixlist.size();
        for (int i = 0; i < nlist; i++) {
            ipix =  pixlist.get(i);
            PixToolsVector3d v = pt.pix2vect_ring(nside,ipix);
            double dist = v.angle(vc);
            assertTrue(dist<=2.*radius);
        }
        cpix = pt.ang2pix_nest(nside,theta,phi);
        PixToolsVector3d vc1 = pt.pix2vect_nest(nside, cpix);
        LongList pixlist1;
        
        radius *=4;
        pixlist1 = new LongList(pt.query_disc_nested(nside, vc1, radius,  inclusive));
        int nlist1 = pixlist1.size();
        for (int i = 0; i < nlist1; i++) {
            ipix =  pixlist1.get(i);
            PixToolsVector3d v = pt.pix2vect_nest(nside,ipix);
            double dist = v.angle(vc1);
            assertTrue(dist<=2.*radius);
        }
        System.out.println(" test query disk  is done -------------------"); 
    }
    /**
     * test Query_disk  check for consistency in the query for RING/NESTED
     */
    public void testQuery_disk2() {
    	System.out.println(" Start test query_disk HiRes!!!!!!!!!!!!!!!!!!!!!!!!");
    	PixToolsNested pt = new PixToolsNested();
    	long nside = 1 << 20 ;
    	double res = pt.PixRes(nside);
    	System.out.println("nside="+nside+" sresolution="+res);
    	double radius = Math.toRadians(res/3600.)/2.;
    	System.out.println("radius="+radius);
    	PixToolsVector3d v1 = new PixToolsVector3d(-0.704, 0.580, 0.408);
    	System.out.println("!!!!!!!!!!!!! NESTED !!!!!!!!!!!");
    	LongRangeSet diskQ = pt.query_disc_nested(nside,
    			v1,
    			radius, true);  // inclusive query at vector point
    	assertEquals("npixels = " + diskQ.size(), 4, diskQ.size() , 1e-1);
    	long pix1 = pt.vect2pix_nest(nside, v1);
    	PixToolsVector3d v2 = pt.pix2vect_nest(nside, pix1);  // vector to pix center
    	//
    	LongRangeSet diskQ2 = pt.query_disc_nested(nside,
    			v2,
    			radius, true);  // inclusive with point at pixel center
    	assertEquals("npixels = " + diskQ2.size(), 5, diskQ2.size() , 1e-1);
    	
    	//
    	LongRangeSet diskQ3 = pt.query_disc_nested(nside,
    			v2,
    			radius, false);  // exclusive with point at pixel center
    	assertEquals("npixels = " + diskQ3.size(), 1, diskQ3.size() , 1e-1);

 //   RING schema   
       	System.out.println("!!!!!!!!!!!!! RING !!!!!!!!!!!");
    	LongRangeSet diskQ4 = pt.query_disc(nside,
    			v1,
    			radius, true);   // inclusiv at vector point 
    	assertEquals("npixels = " + diskQ4.size(), 4, diskQ4.size() , 1e-1);
    	//
  
    	LongRangeSet diskQ5 = pt.query_disc(nside,
    			v2,
    			radius,  true);  // inclusive at pixel center
    	assertEquals("npixels = " + diskQ5.size(), 5, diskQ5.size() , 1e-1);
    	
//    	System.out.println("n pixels in disk5 ="+diskQ5.size());
    	LongRangeSet diskQ6 = pt.query_disc(nside,
    			v2,
    			radius, false);  // exclusive at pixel center
    	assertEquals("npixels = " + diskQ6.size(), 1, diskQ6.size() , 1e-1);
//
//  test HiRes conversions
//    	
    	PixToolsVector3d pos = new PixToolsVector3d( -0.704, 0.580, 0.408 );

    	nside = 1 << 20;
    	System.out.println("HiRes transformation tests: nside="+nside);
    	LongList nestPixels = new LongList(pt.query_disc_nested(nside, pos, radius, true));
    	LongList ringPixels =  new LongList(pt.query_disc(nside, pos, radius, true));
        //sort ringPixels in same order as nestPixels
        for(int i=0;i<ringPixels.size();i++) ringPixels.set(i, pt.ring2nest(nside,ringPixels.get(i)));
        ringPixels = ringPixels.sort();
        for(int i=0;i<ringPixels.size();i++) ringPixels.set(i, pt.nest2ring(nside,ringPixels.get(i)));
    	assertEquals(nestPixels.size(), ringPixels.size());
    	for(int i=0; i< ringPixels.size(); i++) {
    		long iring = ringPixels.get(i);
    		PixToolsVector3d cv = pt.pix2vect_ring(nside, iring);
    		long inest = pt.ring2nest(nside, iring);
    		long inestC = nestPixels.get(i);
    		PixToolsVector3d cvN = pt.pix2vect_nest(nside, inestC);
    		long iringT = pt.nest2ring(nside, inestC);
    		assertEquals(iring,iringT);
    		assertEquals(inest,inestC);
    		assertEquals(" Xv="+cv.x,cv.x,cvN.x,1.e-10);
    		assertEquals(" Yv="+cv.y,cv.y,cvN.y,1.e-10);
    		assertEquals(" Zv="+cv.z,cv.z,cvN.z,1.e-10);
//    		System.out.println(" inest orig="+inestC+" transformed="+inest+" iring orig="+iring+" transf="+iringT);
//    		System.out.println("Vector cv vs cvN x="+cv.x+" cvN.x="+cvN.x);
//    		System.out.println("Vector cv vs cvN y="+cv.y+" cvN.y="+cvN.y);
//    		System.out.println("Vector cv vs cvN z="+cv.z+" cvN.z="+cvN.z);
    		double[] tetphiR = pt.pix2ang_ring(nside, iring);
    		double[] tetphiN= pt.pix2ang_nest(nside, inestC);
    		assertEquals(" theta="+tetphiR[0],tetphiR[0],tetphiN[0],1.e-10);
    		assertEquals(" phi="+tetphiR[1],tetphiR[1],tetphiN[1],1.e-10);
//    		System.out.println("theta R vs N "+tetphiR[0]+" "+tetphiN[0]);
//    		System.out.println("phi R vs N "+tetphiR[1]+" "+tetphiN[1]);
    	}
    	
    	System.out.println(" End test of query_disc2____________________________");
    	
    }
    /**
     * tests GetNside method
     */
    public void testGetNside() {
        System.out.println(" Start test GetNside !!!!!!!!!!!!!!!!!!!!!");

        double pixsize = 0.3;
        PixToolsNested pt = new PixToolsNested();
        long nside = pt.GetNSide(pixsize);
        System.out.println("Requared nside is "+nside);
        assertEquals("nside = " + nside, 1048576, nside, 1e-1);
        System.out.println(" End of GetNSide test _______________________");
    }
    

}

