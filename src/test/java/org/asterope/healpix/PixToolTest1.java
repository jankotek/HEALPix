package org.asterope.healpix;


import junit.framework.TestCase;
import org.apache.commons.math.geometry.Vector3D;


public class PixToolTest1 extends TestCase {

	
	long NSIDE1 = 64;
	
	long NSIDE2 = 512;
	
	public void testCircle1(){
                PixTools tools = new PixTools(NSIDE1);
		Vector3D vect = PixTools.Ang2Vec(0.5178634297507729, 0.06421357206295804);
		System.out.println(vect);
		tools.query_disc( PixTools.Ang2Vec(0.5178634297507729, 0.06421357206295804), 5.817764173314432E-4,  true);
	}
	
	public void testCircle2(){
                PixTools tools = new PixTools(NSIDE1);
		tools.query_disc( PixTools.Ang2Vec(0.3127581538205727, 0.1050979097909252), 0.01454441043328608, true);
	}
	
	public void testCircle3(){
                PixTools tools = new PixTools(NSIDE1);
		Vector3D vect = PixTools.Ang2Vec(0.011620983936195673, 0.44456444930382233);
		System.out.println(vect);
		tools.query_disc( vect, 5.526875964648709E-4, true);
	}

	public void testCircle4(){
                PixTools tools = new PixTools(NSIDE1);
		tools.query_disc( PixTools.Ang2Vec(0.3127581538205727, 0.1050979097909252), 0.01454441043328608,  true);
	}
	
	
	/*** with NSIDE 512 **/	
	public void testCircle5(){
                PixTools tools = new PixTools(NSIDE2);
		tools.query_disc( PixTools.Ang2Vec(1.0486568403767373, 0.036411931519731704), 6.399540590645875E-4,  true);
	}
	   /*** with NSIDE 512 **/ 
    public void testVertexes(){
        double ra = 30.0;
        double dec=30.0;        
        double[] radec = new double[2];
        radec[0] = ra;
        radec[1] = dec;
        double[] thetphi = new double[2];
        thetphi = PixToolsUtils.RaDecToPolar(radec);
        double theta = thetphi[0];
        double phi = thetphi[1];
        PixTools tools = new PixTools(NSIDE2);
        long ipix = tools.ang2pix( theta, phi);
        double[][] vertexes = tools.pix2vertex( ipix);
        long ipixR = tools.ang2pix( theta, phi);
        double[][] vertexesr = tools.pix2vertex( ipixR);
        for (int i=0; i< 4; i++) {
            Vector3D vect = new Vector3D(vertexes[0][i],vertexes[1][i],vertexes[2][i]);
            double[] angs =PixToolsUtils.Vect2Ang(vect);
            Vector3D vectR = new Vector3D(vertexesr[0][i],vertexesr[1][i],vertexesr[2][i]);
            double[] angsR =PixToolsUtils.Vect2Ang(vectR);
            assertEquals("theta="+angs[0],angs[0],angsR[0], 1e-10);
            assertEquals("phi="+angs[1],angs[1],angsR[1], 1e-10);
        }

    }
	   /*** with high res.**/ 
    public void testVertexesHR(){
    	long nside = 1 << 20 ;
        PixTools tools = new PixTools(nside);
        double ra = 30.0;
        double dec=30.0;

        double[] radec = new double[2];
        radec[0] = ra;
        radec[1] = dec;
        double[] thetphi = new double[2];
        thetphi = PixToolsUtils.RaDecToPolar(radec);
        double theta = thetphi[0];
        double phi = thetphi[1];
        long ipix = tools.ang2pix( theta, phi);
        double[][] vertexes = tools.pix2vertex( ipix);
        long ipixR = tools.ang2pix( theta, phi);
        double[][] vertexesr = tools.pix2vertex( ipixR);
        for (int i=0; i< 4; i++) {
            Vector3D vect = new Vector3D(vertexes[0][i],vertexes[1][i],vertexes[2][i]);
            double[] angs =PixToolsUtils.Vect2Ang(vect);
            Vector3D vectR = new Vector3D(vertexesr[0][i],vertexesr[1][i],vertexesr[2][i]);
            double[] angsR =PixToolsUtils.Vect2Ang(vectR);
            assertEquals("theta="+angs[0],angs[0],angsR[0], 1e-10);
            assertEquals("phi="+angs[1],angs[1],angsR[1], 1e-10);
        }
    }
	public void testInverse(){
		for(double ra = 0.0 ; ra < 360.0; ra+=10) {
			for(double de = -85 ; de <= 85; de+=10)
			{	
				
				long area =getHealId(ra,de,NSIDE2);
				RaDePoint p = getRaDe(area,NSIDE2);
		
				//compare with tolerance
				assertTrue(ra+"!="+p.ra,Math.abs(ra- p.ra)<1);
				assertTrue(de+"!="+p.de,Math.abs(de- p.de)<1);
			}		
		}
	}
	
	   /** 
	    * get position ID for ra,de 
	    * @param ra right ascencion in degrees 
	    * @param de declination in degrees
	    * @param nsides - non default number of nsides  
	    */
	   
	    public long getHealId(double ra, double de, long nsides) {
	    	double[] radec = new double[2];
	    	radec[0] = ra;
	    	radec[1] = de;
	    	double[] polar = PixToolsUtils.RaDecToPolar(radec);
                PixTools tools = new PixTools(nsides);
	    	long ip = tools.ang2pix( polar[0], polar[1]);
	    	

	               return ip;
	    }
	
		/** inverse method, convert area ID to RaDePoint */
		public static RaDePoint getRaDe(long area,long nsides) {
                        PixTools tools = new PixTools(nsides);
			double [] polar = tools.pix2ang( area);
			double[]radec = PixToolsUtils.PolarToRaDec(polar);
			return new RaDePoint(radec[0],radec[1]);
		}
		
		
	    /**
	     * Simple point on sky. Used as return value of functions 
	     * <p>
	     * @See CoeliObject for data definition 
	     */
	    static public class RaDePoint {
	        
	    	/** coordinates in degrees*/ 
	    	public double ra,de;
	    	
	            
	    	public RaDePoint() {}
	           
	    	public RaDePoint(double ra, double de){
	    		this.ra = ra;
	    		this.de = de;
	    	}

	            
	    }
}
