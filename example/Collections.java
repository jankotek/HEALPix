import org.asterope.healpix.LongIterator;
import org.asterope.healpix.LongList;
import org.asterope.healpix.LongRangeSet;
import org.asterope.healpix.Region;


/**
 * demonstrates usage of Long collections
 */
public class Collections {
	
	public static void main(String[] args) {
		/* 
		 * This library comes with three collections
		 *  
		 * First is LongRangeSet which is returned by PixTools
		 *   1) is very fast and memory efficient
		 *   2) stores pixel numbers compressed in range
		 *   2) is readonly (and thread safe)
		 *   3) have fast intersect(), union(), complement() operations 
		 *   	which returns new instances of set  
		 *   
		 */
		
		Region r1 = new Region( 20, 20, 30, 30);
		Region r2 = new Region( 30, 30, 30, 30);
		LongRangeSet set1 = r1.pixelize(1024);
		LongRangeSet set2 = r2.pixelize(1024);
		
		//calculate intersection of two areas
		LongRangeSet interset = set1.intersect(set2);
		//calculate union of two areas
		LongRangeSet union = set2.union(set2);
		
		//rest of sky outside of two areas
		LongRangeSet restOfSky = set1.union(set2).complement();
		
		/**
		 * Second collection type is LongSet. Is pretty much as TreeSet<Long>.
		 * Is sorted, modifiable and contains only unique values.
		 * Pixels are not compressed in ranges, but stored in BitSet. 
		 * This is still very effective, because in 1byte it fits 8 pixels   
		 */
		
		//construct new BitSet and fill it with data using iterator
		LongSet b1 = new LongSet();
		LongIterator iter1 = set1.longIterator();
		while(iter1.hasNext()) 
			b1.add(iter1.next());
		
		//or simply use conversion method
		b1 = set1.toLongSet(); 
		
		//add and remove some stuff
		b1.add(100);
		b1.addAll(new long[]{1,2,3,4,5});
		b1.remove(100);
		
		//iterate over all data and sum values
		LongIterator iter2 = b1.longIterator();
		long sum = 0;
		while(iter2.hasNext()){
			sum +=iter2.next();
		}
		
		//and after all operations are done, convert LongSet to more efficient form
		LongRangeSet finalVersion = b1.toLongRangeSet();
	
		/*
		 * third collection type is LongList which corresponds to ArrayList<Long>. 
		 * Is grovable, but more efficient, values are stored in primitive long[] array. 
		 */
		
		LongList l = new LongList();
		l.add(10);
		
	}

}
