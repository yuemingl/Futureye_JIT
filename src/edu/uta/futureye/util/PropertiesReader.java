package edu.uta.futureye.util;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class PropertiesReader {
	Properties p = null;
	
	public PropertiesReader(String fileName) {
		//InputStream inputStream = this.getClass().getClassLoader().getResourceAsStream("./"+fileName);  
		p = new Properties();    
		try {
			InputStream inputStream = new FileInputStream(fileName);    
			p.load(inputStream);
		} catch (IOException e1) {    
			e1.printStackTrace();    
		}
	}
	
	public String getString(String key) {
		return p.getProperty(key);
	}
	
	public Integer getInteger(String key) {
		String v = getString(key);
		if(v == null) return null;
		return Integer.parseInt(v);
	}
	
	public Double getDouble(String key) {
		String v = getString(key);
		if(v == null) return null;
		return Double.parseDouble(v);
	}
	
	public Boolean getBoolean(String key) {
		String v = getString(key);
		if(v == null) return null;
		return Boolean.parseBoolean(v);
	}
	
	public int[] getIntegerArray(String key) {
		String ary = getString(key);
		if(ary == null) return null;
		String[] ss = ary.split(",");
		int[] rlt = new int[ss.length];
		for(int i=0;i<ss.length;i++) {
			rlt[i] = Integer.parseInt(ss[i]);
		}
		return rlt;
	}
	
	public double[] getDoubleArray(String key) {
		String ary = getString(key);
		if(ary == null) return null;
		String[] ss = ary.split(",");
		double[] rlt = new double[ss.length];
		for(int i=0;i<ss.length;i++) {
			rlt[i] = Double.parseDouble(ss[i]);
		}
		return rlt;
	}

//	public static void main(String[] args) {
//		PropertiesReader r = new PropertiesReader("./test.conf");
//		System.out.println(r.getString("A"));
//		System.out.println(r.getString("b"));
//	}

}
