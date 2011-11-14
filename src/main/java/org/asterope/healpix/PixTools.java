//
// Licenced under GPLv2, see licence.txt
// (c) K.M. Gorski, Nickolai Kuropatkin, Jan Kotek,
//


package org.asterope.healpix;


import org.apache.commons.math.geometry.Vector3D;

import java.util.ArrayList;
import java.util.Arrays;


/**
 *
 *  Contains methods translated from HEALPix Fortran90
 *  with increased map resolution in comparison to original Fortran code.
 *  <p>
 *  Is thread safe.
 *
 * @author N Kuropatkin
 * @author Jan Kotek
 *
 */
public class PixTools {


    protected final long nside;

    public PixTools(long nside){
        this.nside = nside;
    }

    public static final class Pixel{
        public final Vector3D centre,north,south,west,east;

        public Pixel(Vector3D centre, Vector3D north, Vector3D south, Vector3D west, Vector3D east) {
            this.centre = centre;
            this.north = north;
            this.south = south;
            this.west = west;
            this.east = east;
        }

        public Vector3D[] toVertex(){
            return new Vector3D[]{west,east,north,south};
        }


    }

    protected static final double twothird = 2. / 3.;

    protected static final double PI = Math.PI;

    protected static final double TWOPI = 2. * PI;

//	private static double FOURPI = 4. * PI;

    protected static final double HALFPI = PI / 2.0;

     static final int ns_max = 1048576; // 2^20


    /**
     * finds pixels having a colatitude (measured from North pole) :
     * theta1 < colatitude < theta2 with 0 <= theta1 < theta2 <= Pi
     * if theta2 < theta1
     * then pixels with 0 <= colatitude < theta2 or theta1 < colatitude < Pi are
     * returned
     *
     * @param theta1
     *            lower edge of the colatitude
     * @param theta2
     *            upper edge of the colatitude
     * @return  LongList of  pixel numbers (long)
     * @throws Exception
     * @throws IllegalArgumentException
     */
    public LongRangeSet query_strip(double theta1, double theta2) throws Exception {
        LongRangeSetBuilder res = new LongRangeSetBuilder();
        long npix, nstrip;
        long iz,  irmin, irmax;
        int is;
        double phi0, dphi;
        double[] colrange = new double[4];
        /* ---------------------------------------- */
        npix = Nside2Npix(nside);
        if (npix < 0) {
            throw new IllegalArgumentException("Nside should be power of 2");
        }
        if ((theta1 < 0.0 || theta1 > PI) || (theta2 < 0.0 || theta2 > PI)) {
            throw new IllegalArgumentException("Illegal value of theta1, theta2");
        }
        if (theta1 <= theta2) {
            nstrip = 1;
            colrange[0] = theta1;
            colrange[1] = theta2;
        } else {
            nstrip = 2;
            colrange[0] = 0.0;
            colrange[1] = theta2;
            colrange[2] = theta1;
            colrange[3] = PI;
        }
        /* loops on strips */
        for (is = 0; is < nstrip; is++) {
            irmin = RingNum(Math.cos(colrange[2 * is]));
            irmax = RingNum(Math.cos(colrange[2 * is + 1]));
            /* loop on ring number */
            for (iz = irmin; iz <= irmax; iz++) {
                phi0 = 0.;
                dphi = PI;
                InRing( iz, phi0, dphi,res);
            }
        }
        return res.build();
    }

    /**
     * finds pixels that lay within a CONVEX polygon defined by its vertex on
     * sphere
     *
     * @param vlist
     *            ArrayList of vectors defining the polygon vertices
     * @param inclusive
     *            if set 1 returns all pixels crossed by polygon boundaries
     * @return  LongList of pixels
     *
     * algorithm: the polygon is divided into triangles vertex 0 belongs to all
     * triangles
     * @throws Exception
     * @throws IllegalArgumentException
     */
    public LongRangeSet query_polygon(ArrayList<Vector3D> vlist,
                                      boolean inclusive) throws Exception {
        LongRangeSet res = LongRangeSetBuilder.EMPTY;
        int nv = vlist.size();
        Vector3D vp0, vp1, vp2;
        Vector3D vo;
//		double surface, fsky;
        double hand;
        double[] ss = new double[nv];
//		int n_in_trg, ilist, ntl;
//        long npix;
        int ix = 0;

        int n_remain, np, nm, nlow;

        //		System.out.println("Start polygon");
        for (int k = 0; k < nv; k++)
            ss[k] = 0.;
        /* -------------------------------------- */
        n_remain = nv;
        if (n_remain < 3) {
            throw new IllegalArgumentException(" Number of vertices should be >= 3");
        }
        /*---------------------------------------------------------------- */
        /* Check that the poligon is convex or has only one concave vertex */
        /*---------------------------------------------------------------- */
        int i0;
        int i2;
        if (n_remain > 3) { // a triangle is always convex
            for (int i1 = 1; i1 <= n_remain - 1; i1++) { // in [0,n_remain-1]
                i0 = (i1 - 1) % n_remain;
                i2 = (i1 + 1) % n_remain;
                vp0 =vlist.get(i0); // select vertices by 3
                // neighbour
                vp1 =vlist.get(i1);
                vp2 = vlist.get(i2);
                // computes handedness (v0 x v2) . v1 for each vertex v1
                vo = Vector3D.crossProduct(vp0,vp2);
                hand = Vector3D.dotProduct(vo,vp1);
                if (hand >= 0.) {
                    ss[i1] = 1.0;
                } else {
                    ss[i1] = -1.0;
                }

            }
            np = 0; // number of vert. with positive handedness
            for (int i = 0; i < nv; i++) {
                if (ss[i] > 0.)
                    np++;
            }
            nm = n_remain - np;

            nlow = Math.min(np, nm);

            if (nlow != 0) {
                if (nlow == 1) { // only one concave vertex
                    if (np == 1) { // ix index of the vertex in the list
                        for (int k = 0; k < nv - 1; k++) {
                            if (Math.abs(ss[k] - 1.0) <= 1.e-12) {
                                ix = k;
                                break;
                            }
                        }
                    } else {
                        for (int k = 0; k < nv - 1; k++) {
                            if (Math.abs(ss[k] + 1.0) <= 1.e-12) {
                                ix = k;
                                break;
                            }
                        }
                    }

                    // rotate pixel list to put that vertex in #0
                    int n_rot = vlist.size() - ix;
                    int ilast = vlist.size() - 1;
                    for (int k = 0; k < n_rot; k++) {
                        Vector3D temp = vlist.get(ilast);
                        vlist.remove(ilast);
                        vlist.add(0,  temp);
                    }
                }
                if (nlow > 1) { // more than 1concave vertex
                    System.out
                            .println(" The polygon has more than one concave vertex");
                    System.out.println(" The result is unpredictable");
                }
            }
        }
        /* fill the polygon, one triangle at a time */
//		npix = Nside2Npix(nside);
        while (n_remain >= 3) {
            vp0 =  vlist.get(0);
            vp1 =  vlist.get(n_remain - 2);
            vp2 =  vlist.get(n_remain - 1);

            /* find pixels within the triangle */
            LongRangeSet templist = query_triangle(vp0, vp1, vp2, inclusive);
            res = res.union(templist);

            n_remain--;
        }

        return res;
    }

    /**
     * generates a list of pixels that lay inside a triangle defined by
     * the three vertex vectors
     *
     * @param v1
     *            Vector3D defines one vertex of the triangle
     * @param v2
     *            Vector3D another vertex
     * @param v3
     *            Vector3D yet another one
     * @param do_inclusive
     *            long 0 (default) only pixels whose centers are inside the
     *            triangle will be listed, if set to 1 all pixels overlaping the
     *            triangle will be listed
     * @return LongList with pixel numbers
     * @throws Exception
     * @throws IllegalArgumentException
     */
    public LongRangeSet query_triangle(
            final Vector3D v1, final Vector3D v2,
            final Vector3D v3,  final boolean  do_inclusive) throws Exception {

        Vector3D[] vv = new Vector3D[3];
        Vector3D[] vo = new Vector3D[3];


        final long npix = Nside2Npix(nside);
        if (npix < 0) {
            throw new IllegalArgumentException(" Nside should be power of 2 >0 and < "+ns_max);
        }
        vv[0] = v1;
        vv[0] = vv[0].normalize();
        vv[1] = v2;
        vv[1] = vv[1].normalize();
        vv[2] = v3;
        vv[2] = vv[2].normalize();


        final double dth1 = 1.0 / (3.0 * (nside * nside));
        final double dth2 = 2.0 / (3.0 * nside);
        /*
                  * determ = (v1 X v2) . v3 determines the left ( <0) or right (>0)
                  * handedness of the triangle
                  */
        Vector3D vt = Vector3D.crossProduct(vv[0],vv[1]);

        final double determ = Vector3D.dotProduct(vt,vv[2]);

        if (Math.abs(determ) < 1.0e-20) {
            throw new IllegalArgumentException("QueryTriangle: the triangle is degenerated - query cannot be performed");
        }
        // The sign of determinant
        final double sdet = determ >= 0. ? 1.0 : -1.0;

        /* vector orthogonal to the great circle containing the vertex doublet */

        vo[0] = Vector3D.crossProduct(vv[1],vv[2]);
        vo[1] = Vector3D.crossProduct(vv[2],vv[0]);
        vo[2] = Vector3D.crossProduct(vv[0],vv[1]);
        vo[0] = vo[0].normalize();
        vo[1] = vo[1].normalize();
        vo[2] = vo[2].normalize();

        /* test presence of poles in the triangle */
        double zmax = -1.0;
        double zmin = 1.0;
        final boolean test1 = (vo[0].getZ() * sdet >= 0.0); // north pole in hemisphere defined by
        // 2-3
        final boolean test2 = (vo[1].getZ() * sdet >= 0.0); // north pole in the hemisphere defined
        // by 1-2
        final boolean test3 = (vo[2].getZ() * sdet >= 0.0); // north pole in hemisphere defined by
        // 1-3
        if (test1 && test2 && test3)
            zmax = 1.0; // north pole in the triangle
        if ((!test1) && (!test2) && (!test3))
            zmin = -1.0; // south pole in the triangle
        /* look for northenest and southernest points in the triangle */


        /* sin of theta for orthogonal vector */
        double[] sto = new double[3];
        for (int i = 0; i < 3; i++) {
            sto[i] = Math.sqrt((1.0 - vo[i].getZ()) * (1.0 + vo[i].getZ()));
        }
        /*
                  * for each segment ( side of the triangle ) the extrema are either -
                  * -the 2 vertices - one of the vertices and a point within the segment
                  */

        zmax = Math.max(Math.max(vv[0].getZ(), vv[1].getZ()), Math.max(vv[2].getZ(), zmax));
        zmin = Math.min(Math.min(vv[0].getZ(), vv[1].getZ()), Math.min(vv[2].getZ(), zmin));
        /*
                  * if we are inclusive, move upper point up, and lower point down, by a
                  * half pixel size
                  */

        double sin_off = 0.0;
        if (do_inclusive) {
            double offset = PI / (4.0 * nside); // half pixel size
            sin_off = Math.sin(offset);
            zmax = Math.min(1.0, Math.cos(Math.acos(zmax) - offset));
            zmin = Math.max(-1.0, Math.cos(Math.acos(zmin) + offset));
        }

        final long irmin = RingNum(zmax);
        final long irmax = RingNum(zmin);

        //		System.out.println("irmin = " + irmin + " irmax =" + irmax);

        /* loop on the rings */
        double[] phi0i = new double[3];
        double[] tgthi = new double[3];
        for (int i = 0; i < 3; i++) {
            tgthi[i] = -1.0e30 * vo[i].getZ();
            phi0i[i] = 0.0;
        }
        for (int j = 0; j < 3; j++) {
            if (sto[j] > 1.0e-10) {
                tgthi[j] = -vo[j].getZ() / sto[j]; // - cotan(theta_orth)

                phi0i[j] = Math.atan2(vo[j].getY(), vo[j].getX()); // Should make it 0-2pi
                // ?
                /* Bring the phi0i to the [0,2pi] domain if need */

                if (phi0i[j] < 0.) {
                    phi0i[j] = BitManipulation.MODULO(
                            (Math.atan2(vo[j].getY(), vo[j].getX()) + TWOPI), TWOPI); //  [0-2pi]
                }

            }
        }
        /*
                  * the triangle boundaries are geodesics: intersection of the sphere
                  * with plans going through (0,0,0) if we are inclusive, the boundaries
                  * are the intersection of the sphere with plains pushed outward by
                  * sin(offset)
                  */
//		double temp = 0.;
        double[][] dom = new double[3][2];
        double[] alldom = new double[6];

        LongRangeSetBuilder res = new LongRangeSetBuilder();
        for (long iz = irmin; iz <= irmax; iz++) {
            boolean found = false;
            final double z;
            if (iz <= nside - 1) { // North polar cap
                z = 1.0 - iz * iz * dth1;
            } else if (iz <= 3 * nside) { // tropical band + equator
                z = (2.0 * nside - iz) * dth2;
            } else {
                z = -1.0 + (4.0 * nside - iz) * (4.0 * nside - iz) * dth1;
            }

            /* computes the 3 intervals described by the 3 great circles */
            final double st = Math.sqrt((1.0 - z) * (1.0 + z));
            final double tgth = z / st; // cotan(theta_ring)
            double[] dc = new double[3];
            for (int j = 0; j < 3; j++) {
                dc[j] = tgthi[j] * tgth - sdet * sin_off
                        / ((sto[j] + 1.0e-30) * st);

            }

            for (int k = 0; k < 3; k++) {

                if (dc[k] * sdet <= -1.0) { // the whole iso-latitude ring is on
                    // right side of the great circle
                    dom[k][0] = 0.0;
                    dom[k][1] = TWOPI;
                } else if (dc[k] * sdet >= 1.0) { // all on the wrong side
                    dom[k][0] = -1.000001 * (k + 1);
                    dom[k][1] = -1.0 * (k + 1);
                } else { // some is good some is bad
                    double phi_neg = phi0i[k] - (Math.acos(dc[k]) * sdet);
                    double phi_pos = phi0i[k] + (Math.acos(dc[k]) * sdet);
                    //
                    if (phi_pos < 0.)
                        phi_pos += TWOPI;
                    if (phi_neg < 0.)
                        phi_neg += TWOPI;

                    //

                    dom[k][0] = BitManipulation.MODULO(phi_neg, TWOPI);
                    dom[k][1] = BitManipulation.MODULO(phi_pos, TWOPI);

                }
                //

            }
            /* identify the intersections (0,1,2 or 3) of the 3 intervals */

            double[] dom12 = PixToolsUtils.intrs_intrv(dom[0], dom[1]);
            final int n12 = dom12.length / 2;
            int ndom = 0;
            if (n12 != 0) {
                if (n12 == 1) {
                    double[]dom123a = PixToolsUtils.intrs_intrv(dom[2], dom12);
                    int n123a = dom123a.length / 2;

                    if (n123a == 0)
                        found = true;
                    if (!found) {
                        System.arraycopy(dom123a, 0, alldom, 0, dom123a.length);

                        ndom = n123a; // 1 or 2
                    }
                }
                if (!found) {
                    if (n12 == 2) {
                        double[] tmp = { dom12[0], dom12[1] };
                        double[]dom123a = PixToolsUtils.intrs_intrv(dom[2], tmp);
                        double[] tmp1 = { dom12[2], dom12[3] };
                        double[]dom123b = PixToolsUtils.intrs_intrv(dom[2], tmp1);
                        int n123a = dom123a.length / 2;
                        int n123b = dom123b.length / 2;
                        ndom = n123a + n123b; // 0, 1, 2 or 3

                        if (ndom == 0)
                            found = true;
                        if (!found) {
                            if (n123a != 0) {
                                System.arraycopy(dom123a, 0, alldom, 0, 2 * n123a);
                            }
                            if (n123b != 0) {
                                for (int l = 0; l < 2 * n123b; l++) {
                                    alldom[l + 2 * n123a] = dom123b[l];
                                }
                            }
                            if (ndom > 3) {
                                throw new InternalError("QueryTriangle: too many intervals found");
                            }
                        }
                    }
                }
                if (!found) {
                    for (long idom = 0; idom < ndom; idom++) {

                        double a_i = alldom[(int) (2 * idom)];
                        double b_i = alldom[(int) (2 * idom + 1)];
                        double phi0 = (a_i + b_i) / 2.0;
                        double dphiring = Math.abs(b_i - a_i) / 2.0;

                        if (dphiring < 0.0) {
                            phi0 += PI;
                            dphiring += PI;
                        }

                        /* finds pixels in the triangle on that ring */
                        InRing(iz, phi0, dphiring,res);
                    }
                }
            }

        }
        return res.build();
    }



    /**
     * generates  all pixels that lays within an
     * angular distance Radius of the center.
     *
     * @param vector
     *            Vector3D pointing to the disc center
     * @param radius
     *            double angular radius of the disc (in RADIAN )
     * @param inclusive
     *            boolean false (default) only pixels whose centers lay in the disc
     *            are listed, if set to true, all pixels overlapping the disc
     *            are listed. In the inclusive mode the radius is increased by half the pixel size.
     *            In this case most probably all neighbor pixels will be listed even with very small
     *            radius.
     *            In case of exclusive search and very small radius when the disc lays completely
     *            inside a pixel the pixel number is returned using vector2pix method.
     * @return  LongList of pixel numbers
     *
     * calls: RingNum(nside, ir) InRing(nside, iz, phi0, dphi,nest) vector2pix(nside,ipix)
     */
    public LongRangeSet query_disc(Vector3D vector, double radius,
                                   boolean inclusive)  {
        LongRangeSetBuilder res = new LongRangeSetBuilder();

        double pixres = PixRes(nside); // in arc seconds
        if (radius < 0.0 || radius > PI) {
            throw new IllegalArgumentException("angular radius is in RADIAN and should be in [0,pi]");
        }

        double dth1 = 1.0 / (3.0 * nside * nside);
        double dth2 = 2.0 / (3.0 * nside);

        double radius_eff = radius;

        //	System.out.println("dth1="+dth1+" dth2="+dth2+" radius="+radius);

        if (inclusive)
            radius_eff += PI / (4.0 * nside); // increase radius by half pixel
        double cosang = Math.cos(radius_eff);
        /* disc center */
        vector = vector.normalize();
        double x0 = vector.getX(); // norm_vect0;
        double y0 = vector.getY(); // norm_vect0;
        double z0 = vector.getZ(); // norm_vect0;
//		System.out.println("x0="+x0+" y0="+y0+" z0="+z0);
        double phi0 = 0.0;
        double dphi = 0.0;
        if (x0 != 0. || y0 != 0.)
            phi0 = BitManipulation.MODULO(Math.atan2(y0, x0) + TWOPI, TWOPI);  // in [0, 2pi]
        double cosphi0 = Math.cos(phi0);
//			System.out.println("phi0="+phi0+" cosphi0="+cosphi0);
        double a = x0 * x0 + y0 * y0;
        /* coordinate z of highest and lowest points in the disc */
        double rlat0 = Math.asin(z0); // latitude in RAD of the center
        double rlat1 = rlat0 + radius_eff;
        double rlat2 = rlat0 - radius_eff;
        double zmax,zmin;
        //
        if (rlat1 >= HALFPI) {
            zmax = 1.0;
        } else {
            zmax = Math.sin(rlat1);
        }

        long irmin = Math.max(1, RingNum(zmax) - 1); // start from a higher point to be safe
        if (rlat2 <= -HALFPI) {
            zmin = -1.0;
        } else {
            zmin = Math.sin(rlat2);
        }

        long irmax = Math.min(4 * nside - 1, RingNum(zmin) + 1); // go down to a lower point
//		System.out.println(" irmax="+irmax+" irmin="+irmin);
//		ilist = -1;
        /* loop on ring number */
        for (long iz = irmin; iz <= irmax; iz++) {
            double z;
            if (iz <= nside - 1) { // north polar cap
                z = 1.0 - iz * iz * dth1;
            } else if (iz <= 3 * nside) { // tropical band + equator
                z = (2.0 * nside - iz) * dth2;
            } else {
                z = -1.0 + (4.0 * nside - iz) * (4.0 * nside - iz) * dth1;
            }
            /* find phi range in the disc for each z */
            double b = cosang - z * z0;
            double c = 1.0 - z * z;
            double cosdphi = b / Math.sqrt(a * c);
            long done = 0;

            if (Math.abs(x0) <= 1.0e-12 && Math.abs(y0) <= 1.0e-12) {
                cosdphi = -1.0;
                dphi = PI;
                done = 1;
            }
            if (done == 0) {
                if (Math.abs(cosdphi) <= 1.0) {
                    dphi = Math.acos(cosdphi); // in [0,pi]
                } else {
                    if (cosphi0 >= cosdphi) {
                        dphi = PI; // all the pixels at this elevation are in
                        // the disc
                    } else {
                        done = 2; // out of the disc
                    }
                }

            }
            if (done < 2) { // pixels in disc
                /* find pixels in the disc */
//				System.out.println("iz="+iz+" phi="+phi0+" dphi="+dphi);

                InRing(iz, phi0, dphi,res);
            }

        }
//
// if no intersections and radius less than pixel size return the pixel number
//
        if (res.size() == 0 && pixres > Math.toDegrees(radius)/3600.) {
            long pixel = vect2pix(vector);

            res.append(pixel);
        }

//make sure that center vector is contained in result set. 
//this is HOTFIX for rounding errors in InRing method
        LongRangeSet ret2 = res.build();
        long centerIpix = vect2pix(vector);
        if(!ret2.contains(centerIpix)){
            //construct new long range set and add it content to result
            LongRangeSetBuilder b2 = new LongRangeSetBuilder();
            b2.append(centerIpix);
            ret2 = ret2.union(b2.build());
        }

        return ret2;


    }

    /**
     * renders theta and phi coordinates of the nominal pixel center for the
     * pixel number ipix (RING scheme) given the map resolution parameter nside
     *
     * @param ipix
     *            long pixel number
     * @return double[] theta,phi
     */
    public double[] pix2ang(long ipix)  {
        double theta,phi;
        /*                            */
        if (nside < 1 || nside > ns_max) {
            throw new IllegalArgumentException("Nside should be power of 2 >0 and < "+ns_max);
        }
        long nsidesq = nside * nside;
        long npix = 12 * nsidesq; // total number of pixels
        if (ipix < 0 || ipix > npix - 1) {
            throw new IllegalArgumentException("ipix out of range calculated from nside");
        }
        long ipix1 = ipix + 1; //  in [1, npix]
        long nl2 = 2 * nside;
        long nl4 = 4 * nside;
        long ncap = 2 * nside * (nside - 1); // points in each polar cap, =0 for
        // nside =1

        if (ipix1 <= ncap) { // North polar cap
            double hip = ipix1 / 2.0;
            double fihip = (long) hip; // get integer part of hip
            long iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from north
            // pole
            long iphi = ipix1 - 2 * iring * (iring - 1);
            theta = Math.acos(1.0 - iring * iring / (3.0 * nsidesq));
            phi = ((double)iphi - 0.5) * PI / (2.0 * iring);


        } else if (ipix1 <= nl2 * (5 * nside + 1)) { // equatorial region
            long ip = ipix1 - ncap - 1;
            long iring = (ip / nl4) + nside; // counted from North pole
            long iphi = ip% nl4 + 1;
            double fodd = 0.5 * (1. + BitManipulation.MODULO(iring + nside, 2)); // 1 if iring+nside
            // is odd, 1/2 otherwise
            theta = Math.acos((nl2 - iring) / (1.5 * nside));
            phi = ((double)iphi - fodd) * PI / (2.0 * nside);

        } else { // South pole cap
            long ip = npix - ipix1 + 1;
            double hip = ip / 2.0;
            double fihip = (long) hip;
            long iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from South
            // pole
            long iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
            theta = Math.acos(-1.0 + iring * iring / (3.0 * nsidesq));
            phi = ((double)iphi - 0.5) * PI / (2.0 * iring);

        }
        double[] res = { 0., 0. };
        res[0] = theta;
        res[1] = phi;
        return res;
    }

    /**
     * returns the vector pointing in the center of the pixel ipix. The vector
     * is calculated by makePix2Vect method
     *
     * @param ipix pixel number
     * @return Vector3D
     */
    public Vector3D pix2vect(long ipix)  {
        return makePix2Vect(ipix).centre;
    }



    /**
     * renders vector (x,y,z) coordinates of the nominal pixel center for pixel
     * ipix (RING scheme) given the map resolution parameter nside. It also
     * calculates (x,y,z) positions of the four vertices in order N,W,S,E. These
     * results are stored in pixVect and pixVertex structures. Those can be
     * obtained using pix2Vect_ring and pix2vert_ring methods
     *
     * @param ipix
     *            pixel number
     * @return Pixel
     */
    public Pixel makePix2Vect(long ipix)  {

        double z_nv, z_sv,  hdelta_phi;
        double   z,phi;
//		boolean do_vertex = true;
        long nsidesq = nside * nside;

        /*                                 */
        if (nside < 1 || nside > ns_max) {
            throw new IllegalArgumentException("Nside should be power of 2 >0 and < "+ns_max);
        }

        long npix = 12 * nsidesq;
        if (ipix < 0 || ipix > npix - 1) {
            throw new IllegalArgumentException("ipix out of range calculated from nside");
        }

        long ipix1 = ipix + 1; //  in [1, npix]
        long nl2 = 2 * nside;
        long nl4 = 4 * nside;
        long ncap = 2 * nside * (nside - 1); // points in each polar cap
        double fact1 = 1.5 * nside;
        double fact2 = 3.0 * nsidesq;
        double phi_nv = 0.0;
        double phi_sv = 0.0;
        if (ipix1 <= ncap) { // north polar cap
            double hip = ipix1 / 2.0;
            double fihip = (long) hip;
            long iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from north
            // pole
            long iphi = ipix1 - 2 * iring * (iring - 1);
            z = 1.0 - iring * iring / fact2;
            phi = (iphi - 0.5) * PI / (2.0 * iring);

            hdelta_phi = PI / (4.0 * iring); // half pixel width
            z_nv = 1.0 - (iring - 1) * (iring - 1) / fact2;
            z_sv = 1.0 - (iring + 1) * (iring + 1) / fact2;
            long iphi_mod =  (iphi - 1) % iring; // in [0,1,...,iring-1]
            long iphi_rat = (iphi - 1) / iring; // in [0,1,2,3]
            if (iring > 1)
                phi_nv = HALFPI * (iphi_rat + iphi_mod / (iring - 1.0));
            phi_sv = HALFPI * (iphi_rat + (iphi_mod + 1.0) / (iring + 1.0));
        } else if (ipix1 <= nl2 * (5 * nside + 1)) { // equatorial region
            long ip =  (ipix1 - ncap - 1);
            long iring = (ip / nl4) + nside; // counted from North pole
            long iphi =  ip% nl4 + 1;
            double fodd = 0.5 * (1. + BitManipulation.MODULO(iring + nside, 2)); // 1 if iring+nside
            // is odd or 1/2
            z = (nl2 - iring) / fact1;
            phi = (iphi - fodd) * PI / (2.0 * nside);
            hdelta_phi = PI / (4.0 * nside); // half pixel width
            phi_nv = phi;
            phi_sv = phi;
            z_nv = (nl2 - iring + 1) / fact1;
            z_sv = (nl2 - iring - 1) / fact1;
            if (iring == nside) { // nothern transition
                z_nv = 1.0 - (nside - 1) * (nside - 1) / fact2;
                long iphi_mod =  (iphi - 1) % nside; // in [0,1,...,nside-1]
                long iphi_rat = (iphi - 1) / nside; // in [0,1,2,3]
                if (nside > 1)
                    phi_nv = HALFPI * (iphi_rat + iphi_mod / (nside - 1.));
            } else if (iring == 3 * nside) { // southern transition
                z_sv = -1.0 + (nside - 1) * (nside - 1) / fact2;
                long iphi_mod = (iphi - 1) % nside; // in [0,1,... iring-1]
                long iphi_rat = (iphi - 1) / nside; // in [0,1,2,3]
                if (nside > 1)
                    phi_sv = HALFPI * (iphi_rat + iphi_mod / (nside - 1.0));
            }

        } else { // South polar cap
            long ip = npix - ipix1 + 1;
            double hip = ip / 2.0;
            double fihip = (long) hip;
            long iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from South
            // pole
            long iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
            z = -1.0 + iring * iring / fact2;
            phi = (iphi - 0.5) * PI / (2.0 * iring);
            hdelta_phi = PI / (4.0 * iring); // half pixel width
            z_nv = -1.0 + (iring + 1) * (iring + 1) / fact2;
            z_sv = -1.0 + (iring - 1) * (iring - 1) / fact2;
            long iphi_mod =  (iphi - 1) % iring; // in [0,1,...,iring-1]
            long iphi_rat = (iphi - 1) / iring; // in [0,1,2,3]
            phi_nv = HALFPI * (iphi_rat + (iphi_mod + 1) / (iring + 1.0));
            if (iring > 1)
                phi_sv = HALFPI * (iphi_rat + iphi_mod / (iring - 1.0));

        }
        /* pixel center */
        double sth = Math.sqrt((1.0 - z) * (1.0 + z));
        //pixVect.getX() = sth * Math.cos(phi);
        //pixVect.getY() = sth * Math.sin(phi);
        //pixVect.getZ() = z;
        Vector3D pixVect = new Vector3D(sth * Math.cos(phi), sth * Math.sin(phi), z);
        /* west vertex */
        double phi_wv = phi - hdelta_phi;

        Vector3D west = new Vector3D(sth * Math.cos(phi_wv), sth * Math.sin(phi_wv),z);

        /* east vertex */
        double phi_ev = phi + hdelta_phi;
        Vector3D east = new Vector3D(sth * Math.cos(phi_ev),  sth * Math.sin(phi_ev), z);
        /* north vertex */
        double sth_nv = Math.sqrt((1.0 - z_nv) * (1.0 + z_nv));
        Vector3D north = new Vector3D(sth_nv * Math.cos(phi_nv), sth_nv * Math.sin(phi_nv), z_nv);
        /* south vertex */
        double sth_sv = Math.sqrt((1.0 - z_sv) * (1.0 + z_sv));
        Vector3D south = new Vector3D(sth_sv * Math.cos(phi_sv), sth_sv * Math.sin(phi_sv),  z_sv);
        return new Pixel(pixVect, north, south, west, east);
    }

    /**
     * renders the pixel number ipix (RING scheme) for a pixel which contains a
     * point with coordinates theta and phi, given the map resolution parameter
     * nside.
     *
     * @param theta
     *            double theta
     * @param phi -
     *            double phi
     * @return  long ipix
     */
    public long ang2pix(double theta, double phi) {

        if (nside < 1 || nside > ns_max) {
            throw new IllegalArgumentException("Nside should be power of 2 >0 and < "+ns_max);
        }
        if (theta < 0.0 || theta > PI) {
            throw new IllegalArgumentException("Theta out of range [0,pi]");
        }

        double z = Math.cos(theta);

        if (phi >= TWOPI)  phi = phi -TWOPI ;

        if (phi < 0.)
            phi =phi + TWOPI; //  phi in [0, 2pi]

        return zPhi2Pix(z, phi);
    }

    /**
     * renders the pixel number ipix (RING scheme) for a pixel which contains a
     * point on a sphere at coordinate vector (x,y,z), given the map resolution
     * parameter nside
     *
     * @param vector
     *            Vector3D of the point coordinates
     * @return  long pixel number
     * @throws IllegalArgumentException
     */
    public long vect2pix(Vector3D vector)  {
        if (nside < 1 || nside > ns_max) {
            throw new IllegalArgumentException("Nside should be power of 2 >0 and < "+ns_max);
        }
        double dnorm = vector.getNorm();
        double z = vector.getZ() / dnorm;
        double phi = 0.;
        if (vector.getX() != 0. || vector.getY() != 0.)
            phi = Math.atan2(vector.getY(), vector.getX()); // phi in [-pi,pi]

        if (phi < 0.)
            phi += TWOPI; //  phi in [0, 2pi]
        return zPhi2Pix(z, phi);
    }

    private long zPhi2Pix(double z, double phi) {
        long ipix1;
        double tt = phi / HALFPI; // in [0,4]
        double za = Math.abs(z);
        long nl2 = 2 * nside;
        long nl4 = 4 * nside;
        long ncap = nl2 * (nside - 1); // number of pixels in the north polar cap
        long npix = 12 * nside * nside;
        if (za < twothird) { // equatorial region
            long jp = (long) (nside * (0.5 + tt - 0.75 * z)); // index of ascending
            // edge line
            long jm = (long) (nside * (0.5 + tt + 0.75 * z)); // index of descending
            // edge line

            long ir = nside + 1 + jp - jm; // in [1,2n+1]
            long kshift = 0;
            if ( BitManipulation.MODULO(ir, 2) == 0)
                kshift = 1; // 1 if ir even, 0 otherwise
            long ip =  ((jp + jm - nside + kshift + 1) / 2) + 1; // in [1,4n]
            ipix1 = ncap + nl4 * (ir - 1) + ip;
        } else { // North and South polar caps
            double tp = tt - (long) tt;
            double tmp = Math.sqrt(3.0 * (1.0 - za));
            long jp = (long) (nside * tp * tmp); // increasing edge line index
            long jm = (long) (nside * (1.0 - tp) * tmp); // decreasing edge index

            long ir = jp + jm + 1; // ring number counted from closest pole
            long ip = (long) (tt * ir) + 1; // in [1,4*ir]
            if (ip > 4 * ir)
                ip = ip - 4 * ir;

            ipix1 = 2 * ir * (ir - 1) + ip;
            if (z <= 0.0)
                ipix1 = npix - 2 * ir * (ir + 1) + ip;
        }
        return ipix1 - 1; // in [0, npix-1]
    }

    public LongRangeSet InRing(long iz, double phi0, double dphi){
        LongRangeSetBuilder b = new LongRangeSetBuilder();
        InRing(iz,phi0,dphi,b);
        return b.build();
    }

    /**
     * returns the list of pixels in RING scheme with latitude in [phi0 -
     * dpi, phi0 + dphi] on the ring iz in [1, 4*nside -1 ] The pixel id numbers
     * are in [0, 12*nside^2 - 1]
     *
     * @param iz
     *           long ring number
     * @param phi0
     *            double
     * @param dphi
     *            double
     * @param res store result in this builder
     *
     * @throws IllegalArgumentException
     *
     * Modified by N. Kuropatkin 07/09/2008  Corrected several bugs and make test of all cases.
     *
     */
    public void InRing(long iz, double phi0, double dphi, LongRangeSetBuilder res)  {

        boolean to_top = false;

        boolean conservative = false;

        final double epsilon = Double.MIN_VALUE; // the constant to eliminate

        long ir = 0;
        final long kshift, nr, ipix1, ipix2,  ncap, npix;//nir1, nir2,
        long ip_low = 0, ip_hi = 0; //,in, nir;
//		long inext;
        npix = 12 * nside * nside; // total number of pixels
        ncap = 2 * nside * (nside - 1); // number of pixels in the north polar
        // cap
        double phi_low = BitManipulation.MODULO(phi0 - dphi, TWOPI) - epsilon; // phi min,
        // excluding
        // 2pi period
        double phi_hi = BitManipulation.MODULO(phi0 + dphi, TWOPI) + epsilon;


        /* identifies ring number */
        if ((iz >= nside) && (iz <= 3 * nside)) { // equatorial region
            ir = iz - nside + 1; // in [1, 2*nside + 1]
            ipix1 = ncap + 4 * nside * (ir - 1); // lowest pixel number in the
            // ring
            ipix2 = ipix1 + 4 * nside - 1; // highest pixel number in the ring
            kshift = (long) BitManipulation.MODULO(ir, 2.);

            nr = nside * 4;
        } else {
            if (iz < nside) { // north pole
                ir = iz;
                ipix1 = 2 * ir * (ir - 1); // lowest pixel number
                ipix2 = ipix1 + 4 * ir - 1; // highest pixel number
            } else { // south pole
                ir = 4 * nside - iz;

                ipix1 = npix - 2 * ir * (ir + 1); // lowest pixel number
                ipix2 = ipix1 + 4 * ir - 1;       // highest pixel number
            }
            nr = ir * 4;
            kshift = 1;
        }

        // Construct the pixel list
        if (Math.abs(dphi - PI) < epsilon) {
            //take entire range
            res.appendRange(ipix1,ipix2);

            return;
        }
        // java calculation jitter
        double shift = kshift / 2.0;

        // conservative : include every intersected pixel, even if the
        // pixel center is out of the [phi_low, phi_hi] region
        if (conservative) {
            ip_low = Math.round((nr * phi_low) / TWOPI - shift);
            ip_hi = Math.round((nr * phi_hi) / TWOPI - shift);

            ip_low = ip_low % nr; // in [0, nr - 1]
            ip_hi = ip_hi% nr; // in [0, nr - 1]
//			System.out.println("ip_low="+ip_low+" ip_hi="+ip_hi);
        } else { // strict: includes only pixels whose center is in
            //                                                    [phi_low,phi_hi]

            ip_low = (long) Math.ceil((nr * phi_low) / TWOPI - shift);
            ip_hi = (long)((nr * phi_hi) / TWOPI - shift);
            if (ip_low == ip_hi + 1)
                ip_low = ip_hi;

            if ((ip_low - ip_hi == 1) && (dphi * nr < PI)) {
                // the interval is too small ( and away from pixel center)
                // so no pixels is included in the list

                System.out
                        .println("the interval is too small and avay from center");

                return; // return empty list
            }

            ip_low = Math.min(ip_low, nr - 1);
            ip_hi = Math.max(ip_hi, 0);
        }

        //
        if (ip_low > ip_hi)
            to_top = true;

        if (to_top) {
            ip_low += ipix1;
            ip_hi += ipix1;
            //nir1 = ipix2 - ip_low + 1;

            //nir2 = ip_hi + 1;
            res.appendRange(ipix1,ip_hi);
            res.appendRange(ip_low,ipix2);
        } else {
            if (ip_low < 0 ){
                ip_low = Math.abs(ip_low) ;
                //nir1 = ip_low;
                //nir2 = ip_hi + 1;

                res.appendRange(ipix1, ipix1+ip_hi);
                res.appendRange(ipix2-ip_low +1, ipix2);

                return ;

            }
            ip_low += ipix1;
            ip_hi += ipix1;

            res.appendRange(ip_low,ip_hi);

        }

    }

    /**
     * returns the ring number in {1, 4*nside - 1} calculated from z coordinate
     *
     * @param z
     *            double z coordinate
     * @return long ring number
     */
    public long RingNum(double z) {
        /* equatorial region */

        long iring = Math.round(nside * (2.0 - 1.5 * z));
        /* north cap */
        if (z > twothird) {
            iring = Math.round(nside * Math.sqrt(3.0 * (1.0 - z)));
            if (iring == 0)
                iring = 1;
        }
        /* south cap */
        if (z < -twothird) {
            iring = Math.round(nside * Math.sqrt(3.0 * (1.0 + z)));
            if (iring == 0)
                iring = 1;
            iring = 4 * nside - iring;
        }
        return iring;
    }

    /**
     * calculates vector corresponding to angles theta (co-latitude
     * measured from North pole, in [0,pi] radians) phi (longitude measured
     * eastward in [0,2pi] radians) North pole is (x,y,z) = (0, 0, 1)
     *
     * @param theta double
     * @param phi double
     * @return Vector3D
     * @throws IllegalArgumentException
     */
    public static Vector3D Ang2Vec(double theta, double phi)  {

        Vector3D v;
        if ((theta < 0.0) || (theta > PI)) {
            throw new IllegalArgumentException("theta out of range [0.,PI]");
        }
        double stheta = Math.sin(theta);
        double x = stheta * Math.cos(phi);
        double y = stheta * Math.sin(phi);
        double z = Math.cos(theta);
        v = new Vector3D(x, y, z);
        return v;
    }

    /**
     * returns nside such that npix = 12*nside^2,  nside should be
     * power of 2 and smaller than ns_max if not return -1
     *
     * @param npix
     *            long the number of pixels in the map
     * @return long nside the map resolution parameter
     */
    public static long Npix2Nside(long npix) {
        long npixmax = 12 *(long) ns_max *(long) ns_max;

        long nside = (long) Math.rint(Math.sqrt(npix / 12));
        if (npix < 12) {
            throw new IllegalArgumentException("npix is too small should be > 12");
        }
        if (npix > npixmax) {
            throw new IllegalArgumentException("npix is too large > 12 * ns_max^2");
        }
        double fnpix = 12.0 * nside * nside;
        if (Math.abs(fnpix - npix) > 1.0e-2) {
            throw new IllegalArgumentException("npix is not 12*nside*nside");
        }
        double flog = Math.log((double) nside) / Math.log(2.0);
        double ilog = Math.rint(flog);
        if (Math.abs(flog - ilog) > 1.0e-6) {
            throw new IllegalArgumentException(" nside is not power of 2");
        }
        return nside;
    }

    /**
     * calculates npix such that npix = 12*nside^2 ,nside should be
     * a power of 2, and smaller than ns_max otherwise return -1
     *
     * @param nside
     *            long the map resolution
     * @return npix long the number of pixels in the map
     */
    public static long Nside2Npix(long nside) {

        long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576,
                2097152, 4194304};

        long res = 0;
        if (Arrays.binarySearch(nsidelist, nside) < 0) {
            throw new IllegalArgumentException("nside should be >0, power of 2, <"+ns_max);
        }
        res = 12 * nside * nside;
        return res;
    }




    /**
     * calculates angular resolution of the pixel map
     * in arc seconds.
     * @param nside
     * @return double resolution in arcsec
     */
    public static double PixRes(long nside) {
        double degrad = Math.toDegrees(1.0);
        double skyArea = 4.*PI*degrad*degrad; // 4PI steredian in deg^2
        double arcSecArea = skyArea*3600.*3600.;  // 4PI steredian in (arcSec^2)
        long npixels = 12*nside*nside;
        double res = arcSecArea/npixels;       // area per pixel
        res = Math.sqrt(res);           // angular size of the pixel arcsec
        return res;
    }
    /**
     * calculate requared nside given pixel size in arcsec
     * @param pixsize in arcsec
     * @return long nside parameter
     */
    public static long GetNSide(double pixsize) {
        long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };
        long res = 0;
        double pixelArea = pixsize*pixsize;
        double degrad = Math.toDegrees(1.);
        double skyArea = 4.*PI*degrad*degrad*3600.*3600.;
        long npixels = (long) (skyArea/pixelArea);
        long nsidesq = npixels/12;
        long nside_req = (long) Math.sqrt(nsidesq);
        long mindiff = ns_max;
        int indmin = 0;
        for (int i=0; i<nsidelist.length; i++) {
            if (Math.abs(nside_req - nsidelist[i]) <= mindiff) {
                mindiff = Math.abs(nside_req - nsidelist[i]);
                res = nsidelist[i];
                indmin = i;
            }
            if ((nside_req > res) && (nside_req < ns_max)) res = nsidelist[indmin+1];
            if (nside_req > ns_max ) {
                System.out.println("nside cannot be bigger than "+ns_max);
                return ns_max;
            }

        }
        return res;
    }

}