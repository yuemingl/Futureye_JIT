package edu.uta.futureye.util;

import org.apache.bcel.generic.ClassGen;

public class FuncClassLoader<T> extends ClassLoader {
	public static FuncClassLoader<?> instance = null;
	
	public FuncClassLoader(ClassLoader classLoader) {
		super(classLoader);
	}
	
	@SuppressWarnings("unchecked")
	public static <T> FuncClassLoader<T> getInstance(ClassLoader classLoader) {
		if(instance == null) {
			FuncClassLoader<T> fcl = new FuncClassLoader<T>(classLoader);
			instance = fcl;
		}
		return (FuncClassLoader<T>) instance;
	}
	
	/**
	 * Return an instance from a ClassGen object
	 *
	 */
	@SuppressWarnings("unchecked")
	public T newInstance(ClassGen cg) {
		byte[] bytes = cg.getJavaClass().getBytes();
		Class<T> cl = null;
		cl = (Class<T>) defineClass(cg.getJavaClass().getClassName(), bytes, 0,
				bytes.length);
		try {
			return cl.newInstance();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Define a class by giving the bytecode array
	 * 
	 * @param name
	 * @param data
	 * @return
	 */
	public Class<?> defineClassForName(String name, byte[] data) {
		return this.defineClass(name, data, 0, data.length);
	}
}
