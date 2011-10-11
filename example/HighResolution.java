import org.asterope.healpix.*;


/** 
 * demonstrates high resolution
 */
public class HighResolution {

	public static void main(String[] args) {
		PixTools t = new PixTools();
		
		/** 
		 * resolution 1 arc second  
		 */
		long nside = t.GetNSide(1);
		System.out.println("NSIDE = "+nside);
		//>> NSIDE = 262144
		
		/**
		 * radius 20 degrees. Result set will have around  1e10 pixels!!
		 */
		double radius = Math.toRadians(20);
		
		/**
		 * random points on equator which is worst case scenario (poles are best case)
		 */		
		PixToolsVector3d p1 =  t.Ang2Vec(Math.toRadians(90), Math.toRadians(30));
		PixToolsVector3d p2 =  t.Ang2Vec(Math.toRadians(90), Math.toRadians(45));
		
		/**
		 * create disc around given pixels
		 */
		LongRangeSet disc1 = t.query_disc(nside, p1, radius, true);
		LongRangeSet disc2 = t.query_disc(nside, p2, radius, true);
		System.out.println("Number of pixels: "+disc1.size()); 
		System.out.println("Number of ranges: "+disc1.rangeCount());
		//>> Number of pixels: 24866171247
		//>> Number of ranges: 268977
		
		
		/**
		 * calculate intersection of two discs
		 */
		LongRangeSet intersect= disc1.intersect(disc2); 
		System.out.println("Intersect number of pixels: "+intersect.size()); 
		System.out.println("Intersect number of ranges: "+intersect.rangeCount());
		//>> Intersect number of pixels: 13374679095
		//>> Intersect number of ranges: 250765

		/*
		 * !!!! and all of this runs in 2 seconds and consumes only 20MB of memory !!!!
		 */
		
	}
}
