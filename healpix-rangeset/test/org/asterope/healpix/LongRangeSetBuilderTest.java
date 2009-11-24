
package org.asterope.healpix;

import junit.framework.TestCase;

public class LongRangeSetBuilderTest extends TestCase{
	
	public void testBuild(){
		LongRangeSetBuilder b = new LongRangeSetBuilder();
		
		b.appendRange(1,5);
		b.appendRange(10,15);
		b.append(16);
		b.appendRange(13,20);
		b.append(21);
		
		LongRangeIterator iter = b.build().rangeIterator();
		assertTrue(iter.moveToNext());
		assertEquals(1,iter.first());
		assertEquals(5,iter.last());

		assertTrue(iter.moveToNext());
		assertEquals(10,iter.first());
		assertEquals(21,iter.last());
		assertTrue(!iter.moveToNext());
	}

}

