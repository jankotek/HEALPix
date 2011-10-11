package org.asterope.healpix;

import junit.framework.TestCase;

import java.util.TreeSet;


public class QueryDiscTest extends TestCase {
	public void testQueryDisc () {
	long nside = 32;
	
	boolean inclusive = false;
	double radius = Math.PI;
	double radius1 = Math.PI/2.;
    PixTools pt = new PixTools(nside);
    int npix = (int) PixTools.Nside2Npix(nside);
    double res = PixTools.PixRes(nside); // pixel size in radians
    System.out.println("res="+res);
    double pixSize = Math.toRadians(res/3600.0); // pixel size in radians
    System.out.println("pixSize="+pixSize+" rad");

    
    LongList fullSky = new LongList(pt.query_disc( new Vector3d(0., 0., 1.), radius,  inclusive));
    LongList firstHalfSky = new LongList(pt.query_disc( new Vector3d(0., 0., 1.), radius1,  inclusive));
    LongList secondHalfSky = new LongList(pt.query_disc( new Vector3d(0., 0., -1.), radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    TreeSet pixHalfsUnique = firstHalfSky.toTreeSet();
    LongList pixHalfsList = new LongList(pixHalfsUnique);
    pixHalfsList = pixHalfsList.sort();
    fullSky = fullSky.sort();

    long listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {

    assertEquals(fullSky.get(i),pixHalfsList.get(i));
    }
    


   firstHalfSky = new LongList(pt.query_disc( new Vector3d(1., 0., 0.), radius1, inclusive));
   secondHalfSky = new LongList(pt.query_disc( new Vector3d(-1., 0., 0.),radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = firstHalfSky.toTreeSet();
    pixHalfsList = new LongList(pixHalfsUnique);
    
    pixHalfsList = pixHalfsList.sort();
    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {
//        System.out.println( "i="+i+" "+fullSky.get(i)+" "+pixHalfsList.get(i));
        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }


    firstHalfSky = new LongList(pt.query_disc( new Vector3d(0., 1., 0.), radius1,  inclusive));
    secondHalfSky = new LongList(pt.query_disc( new Vector3d(0., -1., 0.), radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = firstHalfSky.toTreeSet();
    pixHalfsList = new LongList(pixHalfsUnique);
    pixHalfsList = pixHalfsList.sort();
    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());

    for ( int i=0; i< listL; i++) {

        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }
        
}
}
