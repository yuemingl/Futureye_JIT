package edu.uta.futureye.test.junit;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.VectorEntry;

public class SparseVectorHashMapTest {

	@Test
	public void testIterator() {
		SparseVector v = new SparseVectorHashMap(5);
		int[] idx={1,2,3,5};
		int[] idx2 = new int[4];
		double[] d = {1.0, 2.0, 3.0, 5.0};
		double[] d2 = new double[4];
		v.set(1, d[0]);
		v.set(2, d[1]);
		v.set(3, d[2]);
		v.set(5, d[3]);
		int i=0;
		for(VectorEntry e : v) {
			idx2[i] = e.getIndex();
			d2[i] = e.getValue();
			i++;
		}
		Arrays.sort(idx2);
		Arrays.sort(d2);
//		Assert.assertArrayEquals(d, d2, 0.001);
//		Assert.assertArrayEquals(idx, idx2);
	}

}
