package org.asterope.healpix;


import junit.framework.TestCase;


public class PixToolTest1 extends TestCase {

	static PixToolsNested tools = new PixToolsNested();
	
	long NSIDE1 = 64;
	
	long NSIDE2 = 512;
	
	public void testCircle1(){
		PixToolsVector3d vect = tools.Ang2Vec(0.5178634297507729, 0.06421357206295804);
		System.out.println("vx="+vect.x+" vy="+vect.y+" vz="+vect.z);
		tools.query_disc_nested(NSIDE1, tools.Ang2Vec(0.5178634297507729, 0.06421357206295804), 5.817764173314432E-4,  true);
	}
	
	public void testCircle2(){
		tools.query_disc_nested(NSIDE1, tools.Ang2Vec(0.3127581538205727, 0.1050979097909252), 0.01454441043328608, true);
	}
	
	public void testCircle3(){
		PixToolsVector3d vect = tools.Ang2Vec(0.011620983936195673, 0.44456444930382233);
		System.out.println("vx="+vect.x+" vy="+vect.y+" vz="+vect.z);
		tools.query_disc_nested(NSIDE1, vect, 5.526875964648709E-4, true);
	}

	public void testCircle4(){
		tools.query_disc_nested(NSIDE1, tools.Ang2Vec(0.3127581538205727, 0.1050979097909252), 0.01454441043328608,  true);
	}
	
	
	/*** with NSIDE 512 **/	
	public void testCircle5(){
		tools.query_disc_nested(NSIDE2, tools.Ang2Vec(1.0486568403767373, 0.036411931519731704), 6.399540590645875E-4,  true);
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
        long ipix = tools.ang2pix_nest(NSIDE2, theta, phi);
        double[][] vertexes = tools.pix2vertex_nest(NSIDE2, ipix);
        long ipixR = tools.ang2pix_ring(NSIDE2, theta, phi);
        double[][] vertexesr = tools.pix2vertex_ring(NSIDE2, ipixR);
        for (int i=0; i< 4; i++) {
            PixToolsVector3d vect = new PixToolsVector3d(vertexes[0][i],vertexes[1][i],vertexes[2][i]);
            double[] angs =PixToolsUtils.Vect2Ang(vect);
            PixToolsVector3d vectR = new PixToolsVector3d(vertexesr[0][i],vertexesr[1][i],vertexesr[2][i]);
            double[] angsR =PixToolsUtils.Vect2Ang(vectR);
            assertEquals("theta="+angs[0],angs[0],angsR[0], 1e-10);
            assertEquals("phi="+angs[1],angs[1],angsR[1], 1e-10);
        }

    }
	   /*** with high res.**/ 
    public void testVertexesHR(){
    	long nside = 1 << 20 ;
        double ra = 30.0;
        double dec=30.0;

        double[] radec = new double[2];
        radec[0] = ra;
        radec[1] = dec;
        double[] thetphi = new double[2];
        thetphi = PixToolsUtils.RaDecToPolar(radec);
        double theta = thetphi[0];
        double phi = thetphi[1];
        long ipix = tools.ang2pix_nest(nside, theta, phi);
        double[][] vertexes = tools.pix2vertex_nest(nside, ipix);
        long ipixR = tools.ang2pix_ring(nside, theta, phi);
        double[][] vertexesr = tools.pix2vertex_ring(nside, ipixR);
        for (int i=0; i< 4; i++) {
            PixToolsVector3d vect = new PixToolsVector3d(vertexes[0][i],vertexes[1][i],vertexes[2][i]);
            double[] angs =PixToolsUtils.Vect2Ang(vect);
            PixToolsVector3d vectR = new PixToolsVector3d(vertexesr[0][i],vertexesr[1][i],vertexesr[2][i]);
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
	    	long ip = tools.ang2pix_ring(nsides, polar[0], polar[1]);
	    	

	               return ip;
	    }
	
		/** inverse method, convert area ID to RaDePoint */
		public static RaDePoint getRaDe(long area,long nsides) {
			double [] polar = tools.pix2ang_ring(nsides	, area);
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
