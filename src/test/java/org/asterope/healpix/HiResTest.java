package org.asterope.healpix;

import junit.framework.TestCase;
import org.apache.commons.math.geometry.Vector3D;

//test PixTools on high resolution
public class HiResTest extends TestCase{
	final int NSIDE = 1048576; //highest res available with long ranges
	final double D2R = Math.PI/180d;
	final PixTools ps = new PixTools(NSIDE);
	final Vector3D V = new Vector3D(1, 1, 1).normalize();

	/**
	 * test on high resolutions. If range set are working correctly, no OutOfMemory is generated
	 * 
	 */
	public void testQueryRing() throws Exception{

		
		LongRangeSet rs = null;
		System.out.println("ring 1d");
		rs = ps.query_disc( V, D2R * 1, false);
		System.out.println(rs.size());
		
		System.out.println("ring 10d");
		rs = ps.query_disc( V, D2R * 10, false);
		System.out.println(rs.size());
		
		
	}

}
