package edu.uta.futureye.function;

public enum VN {
	x{ public int getID() {return 0;} },
	y{ public int getID() {return 1;} },
	z{ public int getID() {return 2;} },
	r{ public int getID() {return 3;} },
	s{ public int getID() {return 4;} },
	t{ public int getID() {return 5;} },
	u{ public int getID() {return 6;} },
	v{ public int getID() {return 7;} },
	w{ public int getID() {return 8;} };
	
	public static String[] names = {"x","y","z","r","s","t","u","v","w"};
	
	public abstract int getID();
	
	public static int getID(String ss) {
		if("x".equals(ss)) return 0;
		if("y".equals(ss)) return 1;
		if("z".equals(ss)) return 2;
		if("r".equals(ss)) return 3;
		if("s".equals(ss)) return 4;
		if("t".equals(ss)) return 5;
		if("u".equals(ss)) return 6;
		if("v".equals(ss)) return 7;
		if("w".equals(ss)) return 8;
		return -1;
	}
	
}
