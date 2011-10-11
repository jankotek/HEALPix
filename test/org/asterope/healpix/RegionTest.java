package org.asterope.healpix;

import java.util.ArrayList;

import junit.framework.TestCase;

/**
 * @author kuropat
 * Junit methods to test Region class
 */
public class RegionTest extends TestCase {
	private final static double TWOPI = 2.0 * Math.PI;
//	private final static double epsilon = 1.0e-10;
	/**
	 * test default constructor
	 */
	public void testRegion() {
		double xMin = 0.;
		double xMax = 30.;
		double yMin = 0.;
		double yMax = 10.;
		Region rg = new Region(xMin,xMax,yMin,yMax);
		assertNotNull(rg);
		double x = -10.;
		double y = 15.;
		assertFalse(rg.inReg(x,y));
		x = 10.;
		y = 5.;
		assertTrue(rg.inReg(x,y));
		xMax = -10.;
		ArrayList <Vector3d>vert = rg.getVertices();
		double[][] vertPol = rg.getPolReg();
		for ( int ind=0; ind<vert.size(); ind++) {
			Vector3d vv =  vert.get(ind);
			double [] vvAng = PixToolsUtils.Vect2Ang(vv);
			if (vertPol[ind][1] < 0) vertPol[ind][1] += TWOPI;
//			double comp = BitManipulation.MODULO(vvAng[1], TWOPI) - epsilon;
			assertEquals("theta ",vertPol[ind][0],vvAng[0], 1.0e-5);
			assertEquals("phi ="+vertPol[ind][1],vertPol[ind][1],vvAng[1], 1.0e-5);
		}
		xMin = 20.;
		xMax = 95.;
		Region rg1 = new Region(xMin,xMax,yMin,yMax);
		assertNotNull(rg1);
		x = 45.;
		y = 5.;
		assertTrue(rg1.inReg(x,y));
	}
	/**
	 * test pixelization
	 * 
	 */
	public void testPixelize() {
		System.out.println("test pixelize");
		LongList pixels = new LongList();
		double xMin = 10.;
		double xMax = 60.;
		double yMin = -20.0;
		double yMax = 0.;

		Region rg = new Region(xMin,xMax,yMin,yMax);
		double[][] regCoord = rg.getPolReg();
		for (int i = 0; i<regCoord.length; i++ ) {
			System.out.println("thet="+regCoord[i][0]+" phi="+regCoord[i][1]);
		}
		double resolution = 10.*60.*60.; // resolution in arcsec (= 10 degrees)
		
		try {
			pixels = new LongList(rg.pixelize(resolution));
			long nside = PixTools.GetNSide(resolution);
            PixTools pt = new PixTools(nside);
			int npix = pixels.size();
			assertFalse(npix == 0);
			System.out.println("npix="+npix);
			for (int i=0; i<npix; i++) {
				long pix = ((Long) pixels.get(i)).longValue();
				System.out.println("pixel="+pix);
				double[] pixpos = pt.pix2ang(pix);
				System.out.println("theta="+pixpos[0]+" phi="+pixpos[1]);

				double[][] pixvert = pt.pix2vertex(pix);
				System.out.println("corners");
				for (int j=0; j<pixvert[0].length; j++) {
					double x = pixvert[0][j];
					double y = pixvert[1][j];
					double z = pixvert[2][j];
					double[] pol = new Vector3d(x,y,z).toArray();
					double[] radec1 = PixToolsUtils.PolarToRaDec(pol);
					System.out.println("ra= "+radec1[0]+" dec="+radec1[1]);
				}
				System.out.println();
				
			}
		
		} catch (Exception e) {
			System.err.println("Exception in pixelize");
			e.printStackTrace();
		}

	}
}
