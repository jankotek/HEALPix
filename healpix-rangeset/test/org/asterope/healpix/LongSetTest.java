package org.asterope.healpix;

import junit.framework.TestCase;

public class LongSetTest extends TestCase{
	
	public void testSet(){
		LongSet ls = new LongSet();
		
		assertTrue(ls.isEmpty());
		
		ls.add(100);
		ls.add(101);
		ls.add(102);
		ls.add(Integer.MAX_VALUE);
		
		assertEquals(ls.size(), 4);
		
		LongIterator iter = ls.longIterator();
		assertTrue(iter.hasNext());
		assertEquals(iter.next(), 100);
		assertTrue(iter.hasNext());
		assertEquals(iter.next(), 101);
		assertTrue(iter.hasNext());
		assertEquals(iter.next(), 102);
		assertTrue(iter.hasNext());
		assertEquals(iter.next(),Integer.MAX_VALUE);
		assertTrue(!iter.hasNext());


		
	}

}
