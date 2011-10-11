import org.asterope.healpix.*;


/**
 * Basic example. 
 *
 */
public class HelloWorld {

	public static void main(String[] args) {
		/*
		 * PixTools is central class
		 */
		PixTools t = new PixTools();
		
		/*
		 * decide what arc resolution to use, 
		 * higher resolution means bigger precision, but more calculation and memory usage
		 * maximum precision is 0.4 arc sec. 
		 * use 60 arc sec now
		 */
		long nside = t.GetNSide(60);
		
		/*
		 * points are represented by normalized 3d vectors
		 * get point on sky with RA = 45d, DE=10d. 		  
		 */
		// (80 is because function takes theta instead of DE (distance from north pole))
		double theta = Math.toRadians(80);		
		double phi = Math.toRadians(45);
		Vector3d point = t.Ang2Vec(theta,phi);
		
		//convert point in sky to pixel number
		long ipix = t.vect2pix_ring(nside, point);
		
		//or you can query for pixel number without vector
		ipix = t.ang2pix_ring(nside, theta, phi);
		
		//or query for circular area around point
		//10 degrees radius
		double radius = Math.toRadians(10);
		//include pixels on boundary of circle?
		boolean inclusive = false;
		LongRangeSet area = t.query_disc(nside, point, radius, inclusive);
				
		//check that circle contains original point
		System.out.println("Contains centre: " + area.contains(ipix));
		//check if contains north pole
		long northPoleIpix = t.ang2pix_ring(nside, 0, 0);
		System.out.println("Contains north pole: " + area.contains(northPoleIpix));

	}

}
