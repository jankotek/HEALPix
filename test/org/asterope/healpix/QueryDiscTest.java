package org.asterope.healpix;

import junit.framework.TestCase;


public class QueryDiscTest extends TestCase {
	public void testQueryDisc () {
	long nside = 32;
	
	boolean inclusive = false;
	double radius = Math.PI;
	double radius1 = Math.PI/2.;
    PixTools pt = new PixTools();
    int npix = (int) pt.Nside2Npix(nside);
    double res = pt.PixRes(nside); // pixel size in radians
    System.out.println("res="+res);
    double pixSize = Math.toRadians(res/3600.0); // pixel size in radians
    System.out.println("pixSize="+pixSize+" rad");

    
    LongList fullSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(0., 0., 1.), radius,  inclusive));
    LongList firstHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(0., 0., 1.), radius1,  inclusive));
    LongList secondHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(0., 0., -1.), radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    LongSet pixHalfsUnique = new LongSet(firstHalfSky);
    LongList pixHalfsList = new LongList(pixHalfsUnique);
    pixHalfsList = pixHalfsList.sort();
    fullSky = fullSky.sort();

    long listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {

    assertEquals(fullSky.get(i),pixHalfsList.get(i));
    }
    


   firstHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(1., 0., 0.), radius1, inclusive));
   secondHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(-1., 0., 0.),radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new LongSet(firstHalfSky);
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


    firstHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(0., 1., 0.), radius1,  inclusive));
    secondHalfSky = new LongList(pt.query_disc(nside, new PixToolsVector3d(0., -1., 0.), radius1,  inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new LongSet(firstHalfSky);
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
