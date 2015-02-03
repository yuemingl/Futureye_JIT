package edu.uta.futureye.test.junit;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.VectorEntry;

public class SparseBlockVectorTest {

	@Test
	public void testIterator() {
		//SparseBlockVector 
		SparseBlockVector sbv = new SparseBlockVector(3);
		SparseVector v1 = new SparseVectorHashMap(4);
		SparseVector v2 = new SparseVectorHashMap(1.0,2.0,3.0);
		SparseVector v3 = new SparseVectorHashMap(4.0,5.0);
		sbv.setBlock(3, v3);
		sbv.setBlock(2, v2);
		sbv.setBlock(1, v1);
		sbv.print();
		for(VectorEntry e : sbv) {
			System.out.println(e.getIndex()+" "+e.getValue());
		}

	}

}
