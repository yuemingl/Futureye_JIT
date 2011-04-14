package edu.uta.futureye.test;

import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.ObjVector;

public class ContainerTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ObjList<Integer> list = new ObjList<Integer>();
		list.add(10);
		list.add(20);
		list.add(30);
		list.set(2, 200);
		System.out.println(list);
		
		ObjVector<Integer> vec = new ObjVector<Integer>();
		vec.add(10);
		vec.setSize(5);
		vec.set(3, 30);
		vec.set(5, 50);
		System.out.println(vec);
		
		System.out.println(vec.subVector(2, 4));
		System.out.println(vec.subVector(new ObjIndex(1,3,5)));

	}

}
