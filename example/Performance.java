import org.asterope.healpix.LongRangeSet;
import org.asterope.healpix.PixTools;
import org.asterope.healpix.Vector3d;


/**  
 * measures performance
 */
public class Performance {

	/** helper class to print time */
	static public class StopWatch{
		
		long start = 0;
		public void start(){ 
			start = System.currentTimeMillis();
		}
		
		public void printTime(String label){ 
			long time  = System.currentTimeMillis() - start;
			System.out.println(label+" "+time+" ms");
		}

	}
	
	
	static StopWatch sw = new StopWatch();
	static PixTools t = new PixTools();
	static Vector3d centre = new Vector3d(1,1,1).normalized();
	static long nside = 0;
	static double radius = 0;
	static LongRangeSet result = null;

	public static void main(String[] args) {
						
		nside = t.GetNSide(60);
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
	
		nside = t.GetNSide(60);
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
		
		nside = t.GetNSide(1);
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");

		nside = t.GetNSide(1);
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
		
		nside = 1048576; //highest res available with long ranges
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");

		nside = 1048576; //highest res available with long ranges
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, true);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
	}
}
