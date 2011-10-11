
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

  public void testMoveFirst(){
    LongRangeSetBuilder b = new LongRangeSetBuilder();

    b.appendRange(1,5);
    b.appendRange(10,12);
    b.appendRange(8,16); //builder should extend last range, instead of adding new one

    LongRangeIterator i = b.build().rangeIterator();
    assertTrue(i.moveToNext());
    assertEquals(1, i.first() );
    assertEquals(5, i.last());
    assertTrue(i.moveToNext());
    assertEquals(8, i.first());
    assertEquals(16, i.last());
    assertTrue(!i.moveToNext());

    try{
      //overlap with first range, should throw an exception
      b.appendRange(4,16);
      fail();
    }catch(Exception e){}
  }

}

