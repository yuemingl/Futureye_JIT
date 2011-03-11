package edu.uta.futureye.test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.uta.futureye.util.MultiKey;
import edu.uta.futureye.util.Utils;

public class Test  extends Number implements Comparable<Double> {
    public static Double valueOf(double d) {
        return new Double(d);
    }
	public Test(double d) {
		System.out.println(d);
	}
	
	public static void doubleTest(Test a) {
		System.out.println(a);
	}
	
//	public static void doubleTest(Double a) {
//		System.out.println(a);
//	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//doubleTest(1.0);
		
		// TODO Auto-generated method stub
//		Variable v = new Variable();
//		v.add("x", 3.0);
//		System.out.println(v.getValue("y"));
		
		List<String> a = new LinkedList<String>();
		a.add("x");
		a.add("z");
		List<String> b = new LinkedList<String>();
		b.add("x");
		b.add("y");
		List<String> c= Utils.mergeList(a, b);
		System.out.println(c);
		
//		double eps = 1e-6;
//		double pow = Math.pow(Math.E, -0.25*0.0001/eps);
//		double delta = 0.5*pow/Math.sqrt(Math.PI*eps);
//		System.out.println(delta);
		
//		System.out.println("a".compareTo(null));
		
		Set<Integer> set = new HashSet<Integer>();
		Set<Integer> set2 = new HashSet<Integer>();
		set.add(10);	 
		set.add(20);
		set2.add(20); 
		set2.add(10);
		
		Map<Set,Integer> mapSet = new HashMap<Set,Integer>();
		mapSet.put(set, 10);
		
		System.out.println("-------------");
		System.out.println(set.equals(set2));
		System.out.println(set+" "+set.hashCode());
		System.out.println(set2+" "+set2.hashCode());
		System.out.println(mapSet.get(set2));
		System.out.println("-------------");
	
		
		
		List<Integer> list = new LinkedList<Integer>();
		List<Integer> list2 = new LinkedList<Integer>();
		list.add(1);	 
		list.add(2);
		list2.add(2); 
		list2.add(1);
		
		Map<List,Integer> mapList = new HashMap<List,Integer>();
		mapList.put(list, 100);
		
		System.out.println("-------------");
		System.out.println(list.equals(list2));
		System.out.println(list+" "+list.hashCode());
		System.out.println(list2+" "+list2.hashCode());
		System.out.println(mapList.get(list2));
		System.out.println("-------------");

		
		Map<MultiKey, Integer> map2 = new HashMap<MultiKey, Integer>();
		MultiKey key = new MultiKey(true,1,2);
		map2.put(key,1000);
		System.out.println("-------------");
		System.out.println(map2.get(new MultiKey(true,1,2)));
		System.out.println(map2.get(new MultiKey(true,2,1)));
		System.out.println("-------------");
		
		
		
		
		
	}

	@Override
	public double doubleValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public float floatValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int intValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public long longValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int compareTo(Double o) {
		// TODO Auto-generated method stub
		return 0;
	}

}
