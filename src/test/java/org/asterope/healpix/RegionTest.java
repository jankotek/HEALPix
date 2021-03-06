package org.asterope.healpix;

import junit.framework.TestCase;
import org.apache.commons.math.geometry.Vector3D;

import java.util.ArrayList;

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
		ArrayList <Vector3D>vert = rg.getVertices();
		double[][] vertPol = rg.getPolReg();
		for ( int ind=0; ind<vert.size(); ind++) {
			Vector3D vv =  vert.get(ind);
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

				Vector3D[] pixvert = pt.makePix2Vect(pix).toVertex();
				System.out.println("corners");
				for (int j=0; j<pixvert.length; j++) {
					double x = pixvert[j].getX();
					double y = pixvert[j].getY();
					double z = pixvert[j].getZ();
					double[] pol = new double[]{x,y,z};
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
